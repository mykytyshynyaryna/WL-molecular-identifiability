"""
Benchmark: old (Python-loop) vs new (NumPy-block) flip graph construction.

Usage (from project root):
    python scripts/benchmarks/benchmark_flip_graph.py
    python scripts/benchmarks/benchmark_flip_graph.py --data data/AAAC.smi
    python scripts/benchmarks/benchmark_flip_graph.py --reps 200 --max-mols 500
"""
from __future__ import annotations

import argparse
import itertools
import time
import tracemalloc
from pathlib import Path

import networkx as nx
import numpy as np

from wl_identifiability.flip_graph import build_flip_graph_from_labels



def _group_nodes_by_color(labels: dict) -> dict:
    color2nodes: dict = {}
    for v, c in labels.items():
        color2nodes.setdefault(c, []).append(v)
    return color2nodes



def build_flip_graph(G: nx.Graph, labels: dict) -> nx.Graph:
    A = nx.to_numpy_array(G)
    color_nodes = _group_nodes_by_color(labels)

    for i in color_nodes.values():
        subN = len(i)
        if subN < 2:
            continue
        subA = A[np.ix_(i, i)]
        sub_edgeA = subA.sum()
        if sub_edgeA > (subN * (subN - 1)) / 2:
            flip_matrix = 1 - subA
            np.fill_diagonal(flip_matrix, 0)
            A[np.ix_(i, i)] = flip_matrix

    for i, j in itertools.combinations(color_nodes.values(), 2):
        subB = A[np.ix_(i, j)]
        subN1, subN2 = len(i), len(j)
        if subB.sum() > (subN1 * subN2) / 2:
            flip_matrix = 1 - subB
            A[np.ix_(i, j)] = flip_matrix
            A[np.ix_(j, i)] = flip_matrix.T

    return nx.from_numpy_array(A)



def load_cases(data_path: Path, max_mols: int | None = None) -> list[tuple[str, nx.Graph, dict]]:
    from rdkit import Chem
    from wl_identifiability.graph_construction import convert_rdkit_molecule_to_nx_graph
    from wl_identifiability.wl import compute_wl_coloring

    cases = []
    with data_path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.lower().startswith("smiles"):
                continue
            parts = line.split()
            mol = Chem.MolFromSmiles(parts[0])
            if mol is None:
                continue
            G = convert_rdkit_molecule_to_nx_graph(mol)
            if G.number_of_nodes() == 0:
                continue
            labels = compute_wl_coloring(G, label_attr="atomic_num", store_history=False)["labels"]
            mol_id = parts[1] if len(parts) > 1 else parts[0]
            cases.append((mol_id, G, labels))
            if max_mols is not None and len(cases) >= max_mols:
                break
    return cases



def check_correctness(cases: list[tuple[str, nx.Graph, dict]]) -> bool:
    for label, G, lbls in cases:
        F_numpy = build_flip_graph(G, lbls)
        F_ref, _ = build_flip_graph_from_labels(G, lbls)
        edges_numpy = frozenset(frozenset(e) for e in F_numpy.edges())
        edges_ref   = frozenset(frozenset(e) for e in F_ref.edges())
        if edges_numpy != edges_ref:
            print(f"  CORRECTNESS FAIL on {label}")
            return False
    return True



def _total_time(fn, cases, reps: int) -> float:
    t0 = time.perf_counter()
    for _ in range(reps):
        for _, G, lbls in cases:
            fn(G, lbls)
    return time.perf_counter() - t0


def _peak_memory_kb(fn, cases) -> float:
    tracemalloc.start()
    for _, G, lbls in cases:
        fn(G, lbls)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return peak / 1024



def _build_table(rows: list[dict]) -> str:
    fields = list(rows[0].keys())
    col_w = {f: max(len(f), max(len(str(r[f])) for r in rows)) for f in fields}
    sep    = "+" + "+".join("-" * (col_w[f] + 2) for f in fields) + "+"
    header = "|" + "|".join(f" {f:<{col_w[f]}} " for f in fields) + "|"
    lines  = [sep, header, sep]
    for r in rows:
        lines.append("|" + "|".join(f" {str(r[f]):<{col_w[f]}} " for f in fields) + "|")
    lines.append(sep)
    return "\n".join(lines)



def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data",     default="data/AAAA.smi", help="Path to .smi dataset file")
    parser.add_argument("--max-mols", type=int, default=None,  help="Load at most N molecules")
    parser.add_argument("--reps",     type=int, default=50,    help="Timing repetitions")
    parser.add_argument("--out",      default="profiles",      help="Directory for output table")
    parser.add_argument("--no-save",  action="store_true")
    args = parser.parse_args()

    data_path = Path(args.data)
    if not data_path.exists():
        print(f"ERROR: file not found: {data_path}")
        return

    print(f"Loading molecules from {data_path} ...")
    cases = load_cases(data_path, max_mols=args.max_mols)
    print(f"Loaded {len(cases)} molecules.")

    if not cases:
        print("No valid molecules found.")
        return

    print("Checking correctness ...")
    correct = check_correctness(cases)
    print(f"Correctness: {'PASS' if correct else 'FAIL'}")
    if not correct:
        print("WARNING: timing results may not be meaningful.")

    n_calls = len(cases) * args.reps
    t_numpy = _total_time(build_flip_graph,          cases, args.reps)
    t_ref   = _total_time(build_flip_graph_from_labels, cases, args.reps)
    mem_numpy = _peak_memory_kb(build_flip_graph,          cases)
    mem_ref   = _peak_memory_kb(build_flip_graph_from_labels, cases)

    speedup_time = t_ref / t_numpy if t_numpy > 0 else float("inf")
    speedup_mem  = mem_ref / mem_numpy if mem_numpy > 0 else float("inf")

    rows = [
        {
            "implementation": "ref (loops)",
            "total_s":        f"{t_ref:.4f}",
            "ms/call":        f"{1000 * t_ref / n_calls:.4f}",
            "peak_mem_KB":    f"{mem_ref:.1f}",
        },
        {
            "implementation": "new (numpy)",
            "total_s":        f"{t_numpy:.4f}",
            "ms/call":        f"{1000 * t_numpy / n_calls:.4f}",
            "peak_mem_KB":    f"{mem_numpy:.1f}",
        },
        {
            "implementation": "speedup",
            "total_s":        f"{speedup_time:.2f}x",
            "ms/call":        f"{speedup_time:.2f}x",
            "peak_mem_KB":    f"{speedup_mem:.2f}x",
        },
    ]

    table = _build_table(rows)
    print(f"\nDataset: {data_path}   Molecules: {len(cases)}   Reps: {args.reps}   Total calls: {n_calls}\n")
    print(table)

    if not args.no_save:
        out_dir = Path(args.out)
        out_dir.mkdir(parents=True, exist_ok=True)
        stem = f"flip_graph_benchmark_{data_path.stem}"
        out_path = out_dir / f"{stem}.txt"
        out_path.write_text(table + "\n", encoding="utf-8")
        print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
