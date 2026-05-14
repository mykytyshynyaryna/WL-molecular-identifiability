from __future__ import annotations

import sys
from multiprocessing import Pool
from typing import Any

import networkx as nx
import numpy as np
import polars as pl
from rdkit import Chem

from .bouquet import analyze_bouquet_forest_structure
from .flip_graph import build_flip_graph_from_labels
from .graph_construction import convert_rdkit_molecule_to_nx_graph
from .skeleton import REL_IRREGULAR, build_skeleton, check_lemma16_conditions
from .wl import compute_fixed_wl_steps_from_topology, compute_wl_with_fixed_steps


def estimate_global_wl_steps(graphs: list[nx.Graph], cap: int = 50) -> dict[str, Any]:
    """
    Estimate a global number of WL iterations across a dataset.
    """
    conv_iters = []

    for G in graphs:
        K = compute_fixed_wl_steps_from_topology(G, cap=cap)
        conv_iters.append(K)

    K_max = int(max(conv_iters)) if conv_iters else 0
    K_p95 = int(np.percentile(conv_iters, 95)) if conv_iters else 0

    return {
        "K_max": K_max,
        "K_p95": K_p95,
        "all_iters": conv_iters,
    }


def _collect_sample_graphs_for_k_estimation(
    df: pl.DataFrame,
    sample_size: int = 300,
) -> list[nx.Graph]:
    """Parse the first sample_size SMILES strings from df and return them as NetworkX graphs."""
    sample_graphs: list[nx.Graph] = []

    for smiles in df["smiles"][:sample_size].to_list():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        graph = convert_rdkit_molecule_to_nx_graph(mol)
        sample_graphs.append(graph)

    return sample_graphs


def estimate_fixed_wl_steps_from_dataframe(
    df: pl.DataFrame,
    sample_size: int = 300,
    cap: int = 50,
) -> dict[str, Any]:
    """
    Estimate a fixed number of WL iterations from a sample of molecular graphs.
    """
    sample_graphs = _collect_sample_graphs_for_k_estimation(df, sample_size=sample_size)
    return estimate_global_wl_steps(sample_graphs, cap=cap)


def _compute_wl_summary(wl_result: dict[str, Any], prefix: str) -> dict[str, Any]:
    """Extract convergence, iteration count, and number of color classes from a WL result dict.

    Fields saved:
      wl_converged_{prefix} — True if partition stabilised before the budget was exhausted.
      wl_iters_{prefix}     — Step at which convergence was detected (the actual refinement
                              depth), or None if the full budget was exhausted without
                              convergence.  This is the semantically meaningful value.
      wl_budget_{prefix}    — Total WL steps executed (= K in fixed-step mode).  Every
                              molecule in a run shares the same budget; stored here for
                              per-row traceability.
      n_colors_{prefix}     — Number of distinct color classes in the final labelling.
    """
    labels = wl_result.get("labels", {})
    n_colors = len(set(labels.values())) if labels else 0

    return {
        f"wl_converged_{prefix}": wl_result.get("converged"),
        f"wl_iters_{prefix}": wl_result.get("converge_iter"),
        f"wl_budget_{prefix}": wl_result.get("iterations"),
        f"n_colors_{prefix}": n_colors,
    }


def _compute_flip_graph_summary(flip_info: dict[str, Any]) -> dict[str, Any]:
    """Flatten the flip_info counters into a flat dict suitable for a result row."""
    within_copy = flip_info.get("within_copy", 0) or 0
    within_flip = flip_info.get("within_flip", 0) or 0
    between_copy = flip_info.get("between_copy", 0) or 0
    between_flip = flip_info.get("between_flip", 0) or 0

    return {
        "within_copy": within_copy,
        "within_flip": within_flip,
        "between_copy": between_copy,
        "between_flip": between_flip,
        "n_flipped_edges": within_flip + between_flip,
    }


def _compute_color_ratio(wl_top: dict[str, Any], wl_atom: dict[str, Any]) -> float | None:
    """Return the ratio of atom-aware color classes to topological color classes, or None if undefined."""
    top_labels = wl_top.get("labels", {})
    atom_labels = wl_atom.get("labels", {})

    n_top = len(set(top_labels.values())) if top_labels else 0
    n_atom = len(set(atom_labels.values())) if atom_labels else 0

    if n_top == 0:
        return None

    return n_atom / n_top


def _compute_skeleton_summary(
    class_relations: dict[Any, Any],
    within_class_structures: dict[Any, str],
) -> dict[str, Any]:
    """
    Summarise the skeleton for the result dict.

    Counts Lemma 13 violations (within-class 'other' structures) and
    Lemma 14 violations (irregular between-class relations) so the caller
    can quickly identify theory violations without parsing the full skeleton.
    """
    n_within_other = sum(1 for s in within_class_structures.values() if s == "other")
    n_irregular = sum(1 for rel, _k, _l in class_relations.values() if rel == REL_IRREGULAR)
    return {
        "skeleton_n_classes": len(within_class_structures),
        "n_within_other": n_within_other,
        "n_irregular_class_pairs": n_irregular,
    }


def analyze_single_molecule(smiles: str, molecule_id: str, fixed_wl_steps: int) -> dict[str, Any]:
    """
    Run the full identifiability analysis pipeline for a single molecule.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"[WARN] parse failed: molecule_id={molecule_id} smiles={smiles}", file=sys.stderr)
        return {
            "molecule_id": molecule_id,
            "smiles": smiles,
            "ok": False,
            "stage": "parse",
            "reason": "rdkit_parse_failed",
        }

    graph = convert_rdkit_molecule_to_nx_graph(mol)
    n_nodes = graph.number_of_nodes()
    n_edges = graph.number_of_edges()

    try:
        wl_top = compute_wl_with_fixed_steps(graph, fixed_wl_steps, label_attr=None)
        wl_atom = compute_wl_with_fixed_steps(graph, fixed_wl_steps, label_attr="atomic_num")

    except Exception as error:
        print(
            f"[WARN] wl failed: molecule_id={molecule_id}  {type(error).__name__}: {error}",
            file=sys.stderr,
        )
        return {
            "molecule_id": molecule_id,
            "smiles": smiles,
            "ok": False,
            "stage": "wl",
            "reason": f"wl_error: {type(error).__name__}: {error}",
            "n_nodes": n_nodes,
            "n_edges": n_edges,
        }

    labels_top = wl_top["labels"]
    labels_atom = wl_atom["labels"]

    wl_top_summary = _compute_wl_summary(wl_top, "top")
    wl_atom_summary = _compute_wl_summary(wl_atom, "atom")

    try:
        flip_graph_atom, flip_info_atom = build_flip_graph_from_labels(graph, labels_atom)
        flip_graph_top, _flip_info_top = build_flip_graph_from_labels(graph, labels_top)

    except Exception as error:
        print(
            f"[WARN] flip_graph failed: molecule_id={molecule_id}  {type(error).__name__}: {error}",
            file=sys.stderr,
        )
        return {
            "molecule_id": molecule_id,
            "smiles": smiles,
            "ok": False,
            "stage": "flip_graph",
            "reason": f"flip_error: {type(error).__name__}: {error}",
            "n_nodes": n_nodes,
            "n_edges": n_edges,
            **wl_top_summary,
            **wl_atom_summary,
        }

    flip_summary = _compute_flip_graph_summary(flip_info_atom)

    try:
        skeleton, class_relations, within_structures = build_skeleton(flip_graph_atom, labels_atom)
        lemma16 = check_lemma16_conditions(skeleton, within_structures)
        skeleton_summary = _compute_skeleton_summary(class_relations, within_structures)

    except Exception as error:
        print(
            f"[WARN] skeleton failed: molecule_id={molecule_id}  {type(error).__name__}: {error}",
            file=sys.stderr,
        )
        return {
            "molecule_id": molecule_id,
            "smiles": smiles,
            "ok": False,
            "stage": "skeleton",
            "reason": f"skeleton_error: {type(error).__name__}: {error}",
            "n_nodes": n_nodes,
            "n_edges": n_edges,
            **wl_top_summary,
            **wl_atom_summary,
            **flip_summary,
        }

    def _run_bouquet(flip_graph: nx.Graph, labels: dict[Any, Any], mode_tag: str) -> dict[str, Any]:
        try:
            return analyze_bouquet_forest_structure(flip_graph, labels=labels)
        except Exception as error:
            print(
                f"[WARN] bouquet_forest({mode_tag}) failed: molecule_id={molecule_id}  {type(error).__name__}: {error}",
                file=sys.stderr,
            )
            return {
                "is_bouquet_forest": None,
                "reason": "unknown_error",
            }

    bf_top = _run_bouquet(flip_graph_top, labels_top, "top")
    bf_atom = _run_bouquet(flip_graph_atom, labels_atom, "atom")

    def _verdict(bf: dict[str, Any]) -> int | None:
        v = bf.get("is_bouquet_forest")
        return int(v) if v is not None else None

    return {
        "molecule_id": molecule_id,
        "smiles": smiles,
        "ok": True,
        "stage": "done",
        "n_nodes": n_nodes,
        "n_edges": n_edges,
        **wl_top_summary,
        **wl_atom_summary,
        "color_ratio_atom_to_top": _compute_color_ratio(wl_top, wl_atom),
        **flip_summary,
        **skeleton_summary,
        "lemma16_all_satisfied": lemma16.get("all_satisfied"),
        "lemma16_violations": lemma16.get("violations", []),
        "top_bouquet_forest_verdict": _verdict(bf_top),
        "top_bouquet_forest_reason": bf_top.get("reason"),
        "top_has_bouquet_component": bf_top.get("has_bouquet_component"),
        "atom_bouquet_forest_verdict": _verdict(bf_atom),
        "atom_bouquet_forest_reason": bf_atom.get("reason"),
        "atom_has_bouquet_component": bf_atom.get("has_bouquet_component"),
    }


def _worker_task(args: tuple[str, str, int]) -> dict[str, Any]:
    """Unpacks args and calls analyze_single_molecule. Must be top-level for pickling."""
    smiles, molecule_id, fixed_wl_steps = args
    return analyze_single_molecule(smiles=smiles, molecule_id=molecule_id, fixed_wl_steps=fixed_wl_steps)


def run_molecule_analysis_pipeline(
    df: pl.DataFrame,
    fixed_wl_steps: int,
    n_workers: int = 1,
) -> tuple[pl.DataFrame, int]:
    """
    Run the full molecule analysis pipeline for all rows in the dataframe.

    Returns a Polars DataFrame of results and the count of failed molecules.
    SQLite persistence is the responsibility of the calling script.
    """
    tasks = [
        (smiles, molecule_id, fixed_wl_steps)
        for smiles, molecule_id in zip(df["smiles"].to_list(), df["molecule_id"].to_list(), strict=True)
    ]

    rows: list[dict[str, Any]] = []
    bad = 0

    if n_workers > 1:
        with Pool(processes=n_workers) as pool:
            for row in pool.imap_unordered(_worker_task, tasks, chunksize=1):
                rows.append(row)
                if not row["ok"]:
                    bad += 1
    else:
        for task in tasks:
            row = _worker_task(task)
            rows.append(row)
            if not row["ok"]:
                bad += 1

    result_df = pl.DataFrame(rows)
    return result_df, bad
