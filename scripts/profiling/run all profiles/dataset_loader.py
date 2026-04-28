"""
Shared dataset loader for profiling scripts.

Reads a ZINC-format .smi file (space-separated, with 'smiles zinc_id' header),
parses each SMILES into a NetworkX graph, and runs 1-WL to produce node labels.

Functions
---------
load_cases(data_path, max_mols)
    Load molecules from a single .smi file.

load_cases_from_dir(data_dir, max_mols)
    Load and combine molecules from all .smi files in a directory.
    Returns (cases, n_files).
"""
from __future__ import annotations

from pathlib import Path

import networkx as nx


def load_cases(
    data_path: Path | str,
    max_mols: int | None = None,
) -> list[tuple[str, nx.Graph, dict]]:
    """
    Load molecules from a .smi file and return profiling cases.

    Each case is a (mol_id, G, labels) triple where:
      - mol_id  : ZINC identifier string
      - G       : NetworkX graph with 'atomic_num' node attributes
      - labels  : dict[node -> wl_color] after 1-WL convergence

    Parameters
    ----------
    data_path : path to a .smi file (space-separated, smiles zinc_id columns)
    max_mols  : stop after this many successfully parsed molecules (None = all)
    """
    from rdkit import Chem
    from wl_identifiability.graph_construction import convert_rdkit_molecule_to_nx_graph
    from wl_identifiability.wl import compute_wl_coloring

    data_path = Path(data_path)
    cases: list[tuple[str, nx.Graph, dict]] = []

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
            labels = compute_wl_coloring(
                G, label_attr="atomic_num", store_history=False
            )["labels"]
            mol_id = parts[1] if len(parts) > 1 else parts[0]
            cases.append((mol_id, G, labels))
            if max_mols is not None and len(cases) >= max_mols:
                break

    return cases


def load_cases_from_dir(
    data_dir: Path | str,
    max_mols: int | None = None,
) -> tuple[list[tuple[str, nx.Graph, dict]], int]:
    """
    Load and combine molecules from all .smi files in a directory.

    Files are processed in sorted order. Loading stops once max_mols total
    molecules have been collected (across all files).

    Returns
    -------
    (cases, n_files) where:
      - cases   : combined list of (mol_id, G, labels) triples
      - n_files : number of .smi files that were read
    """
    data_dir = Path(data_dir)
    smi_files = sorted(data_dir.glob("*.smi"))
    if not smi_files:
        return [], 0

    all_cases: list[tuple[str, nx.Graph, dict]] = []
    n_files = 0

    for smi_path in smi_files:
        remaining = None if max_mols is None else max_mols - len(all_cases)
        if remaining is not None and remaining <= 0:
            break
        batch = load_cases(smi_path, max_mols=remaining)
        all_cases.extend(batch)
        n_files += 1

    return all_cases, n_files
