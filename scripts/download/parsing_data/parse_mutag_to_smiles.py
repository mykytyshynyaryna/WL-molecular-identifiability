"""
Parse the MUTAG dataset (TU-format) and convert each graph to a SMILES string.

Node labels (atom types):
  0 -> C, 1 -> N, 2 -> O, 3 -> F, 4 -> I, 5 -> Cl, 6 -> Br

Edge labels (bond types):
  0 -> aromatic, 1 -> single, 2 -> double, 3 -> triple

Output: data/processed/MUTAG/mutag_smiles.smi
  A plain-text, space-separated file with exactly two columns:
      smiles zinc_id
  where zinc_id is the graph index (1-based) from the TU dataset.
  Rows with invalid or missing SMILES are silently skipped.
  Encoding: UTF-8, Unix line endings.

How to run (from the project root):
    python scripts/download/parse_mutag_to_smiles.py

Or with explicit paths:
    python scripts/download/parse_mutag_to_smiles.py \\
        --data-dir data/raw/MUTAG/MUTAG \\
        --out data/processed/MUTAG/mutag_smiles.smi
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from collections import defaultdict

ROOT_DIR = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT_DIR))

from rdkit import Chem
from rdkit.Chem import RWMol, SanitizeMol

ATOM_MAP: dict[int, str] = {0: "C", 1: "N", 2: "O", 3: "F", 4: "I", 5: "Cl", 6: "Br"}

BOND_MAP: dict[int, Chem.rdchem.BondType] = {
    0: Chem.rdchem.BondType.AROMATIC,
    1: Chem.rdchem.BondType.SINGLE,
    2: Chem.rdchem.BondType.DOUBLE,
    3: Chem.rdchem.BondType.TRIPLE,
}



def _read_column(path: Path) -> list[int]:
    """Read a single-column integer file (one value per line)."""
    return [int(line.strip()) for line in path.read_text().splitlines() if line.strip()]


def _read_edges(path: Path) -> list[tuple[int, int]]:
    """Read adjacency list: each line is 'u, v' (1-indexed node ids)."""
    edges = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        u, v = [int(x) for x in line.split(",")]
        edges.append((u, v))
    return edges


def parse_mutag(data_dir: Path) -> list[dict]:
    """
    Parse MUTAG TU-format files and return a list of graph dicts:
        {graph_id, label, nodes: {node_id: atom_symbol}, edges: [(u, v, bond_type_int)]}
    """
    graph_indicator = _read_column(data_dir / "MUTAG_graph_indicator.txt")
    node_labels     = _read_column(data_dir / "MUTAG_node_labels.txt")
    graph_labels    = _read_column(data_dir / "MUTAG_graph_labels.txt")
    edges           = _read_edges(data_dir / "MUTAG_A.txt")
    edge_labels     = _read_column(data_dir / "MUTAG_edge_labels.txt")

    nodes_per_graph: dict[int, dict[int, str]] = defaultdict(dict)
    for node_id, (graph_id, label_idx) in enumerate(
        zip(graph_indicator, node_labels), start=1
    ):
        nodes_per_graph[graph_id][node_id] = ATOM_MAP[label_idx]

    edges_per_graph: dict[int, list[tuple[int, int, int]]] = defaultdict(list)
    for (u, v), bond_label in zip(edges, edge_labels):
        if u < v:
            graph_id = graph_indicator[u - 1]
            edges_per_graph[graph_id].append((u, v, bond_label))

    graphs = []
    for graph_id, label in enumerate(graph_labels, start=1):
        graphs.append(
            {
                "graph_id": graph_id,
                "label": label,
                "nodes": nodes_per_graph[graph_id],
                "edges": edges_per_graph[graph_id],
            }
        )
    return graphs



_BOND_ORDER = {
    Chem.rdchem.BondType.SINGLE:   1,
    Chem.rdchem.BondType.DOUBLE:   2,
    Chem.rdchem.BondType.TRIPLE:   3,
    Chem.rdchem.BondType.AROMATIC: 1,
}

_NORMAL_VALENCE: dict[str, int] = {"C": 4, "N": 3, "O": 2, "F": 1,
                                    "Cl": 1, "Br": 1, "I": 1}


def _explicit_valence(atom: Chem.rdchem.Atom) -> int:
    """Sum of bond orders for all bonds of this atom (pre-sanitization)."""
    return sum(
        _BOND_ORDER.get(bond.GetBondType(), 1)
        for bond in atom.GetBonds()
    )


def _assign_formal_charges(mol: RWMol) -> None:
    """
    Heuristically fix formal charges for atoms that exceed their normal valence.

    MUTAG encodes nitro groups (-NO2) without formal charges: the N has
    explicit bond-order sum of 4 (single to ring C + double to O + single to O)
    while RDKit allows neutral N a maximum of 3.  We detect this and apply the
    standard nitro representation: N(+1) / O(-1).
    """
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N":
            continue
        if atom.GetFormalCharge() != 0:
            continue
        ev = _explicit_valence(atom)
        normal = _NORMAL_VALENCE.get("N", 3)
        if ev > normal:
            atom.SetFormalCharge(1)
            o_minus_assigned = False
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                if (
                    neighbor.GetSymbol() == "O"
                    and neighbor.GetDegree() == 1
                    and bond.GetBondType() == Chem.rdchem.BondType.SINGLE
                    and not o_minus_assigned
                ):
                    neighbor.SetFormalCharge(-1)
                    o_minus_assigned = True


def graph_to_smiles(graph: dict) -> str | None:
    """
    Build an RDKit RWMol from node/edge data and return a canonical SMILES string.
    Returns None if sanitization fails.
    """
    nodes: dict[int, str] = graph["nodes"]
    edges: list[tuple[int, int, int]] = graph["edges"]

    sorted_node_ids = sorted(nodes.keys())
    node_to_idx: dict[int, int] = {nid: i for i, nid in enumerate(sorted_node_ids)}

    mol = RWMol()

    for nid in sorted_node_ids:
        atom = Chem.Atom(nodes[nid])
        mol.AddAtom(atom)

    for u, v, bond_label in edges:
        bond_type = BOND_MAP[bond_label]
        mol.AddBond(node_to_idx[u], node_to_idx[v], bond_type)

    _assign_formal_charges(mol)

    try:
        SanitizeMol(mol)
        return Chem.MolToSmiles(mol)
    except Exception:
        try:
            SanitizeMol(
                mol,
                Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES,
            )
            smiles = Chem.MolToSmiles(mol)
            if Chem.MolFromSmiles(smiles) is None:
                print(
                    f"  [warn] graph_id={graph['graph_id']} relaxed SMILES failed "
                    f"round-trip (over-valent atom at fused-ring junction): {smiles}"
                )
                return None
            print(f"  [info] graph_id={graph['graph_id']} used relaxed sanitization: {smiles}")
            return smiles
        except Exception as exc2:
            print(f"  [warn] graph_id={graph['graph_id']} sanitization failed: {exc2}")
            return None



def main() -> None:
    parser = argparse.ArgumentParser(description="Convert MUTAG dataset to SMILES (.smi).")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=ROOT_DIR / "data" / "raw" / "MUTAG" / "MUTAG",
        help="Directory containing MUTAG_*.txt files",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT_DIR / "data" / "processed" / "MUTAG" / "mutag_smiles.smi",
        help="Output .smi path (space-separated: smiles zinc_id)",
    )
    args = parser.parse_args()

    data_dir: Path = args.data_dir
    out_path: Path = args.out

    print(f"Parsing MUTAG from: {data_dir}")
    graphs = parse_mutag(data_dir)
    print(f"Found {len(graphs)} graphs")

    rows: list[tuple[str, int]] = []
    failed = 0
    for g in graphs:
        smiles = graph_to_smiles(g)
        if smiles is None:
            failed += 1
        else:
            rows.append((smiles, g["graph_id"]))

    print(f"Converted {len(rows)}/{len(graphs)} graphs successfully ({failed} failed)")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="\n") as f:
        f.write("smiles zinc_id\n")
        for smiles, zinc_id in rows:
            f.write(f"{smiles} {zinc_id}\n")

    print(f"Saved to: {out_path}")

    print("\nSample (first 5 rows):")
    print("  smiles zinc_id")
    for smiles, zinc_id in rows[:5]:
        print(f"  {smiles} {zinc_id}")


if __name__ == "__main__":
    main()
