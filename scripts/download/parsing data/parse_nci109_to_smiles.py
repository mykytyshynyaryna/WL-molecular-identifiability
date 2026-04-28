"""
Parse the NCI109 dataset (TU-format) and convert each graph to a SMILES string.

NCI109 differences from MUTAG
------------------------------
* No edge_labels.txt — bond types are not stored.  All bonds are built as
  SINGLE and valence is satisfied through implicit hydrogens and (where
  necessary) formal charges on over-valent atoms.
* Node labels are 1-indexed integers (1–38).  NCI109 is the companion dataset
  to NCI1 and uses the same encoding scheme (labels 1–37) plus label 38 = Se:

    label 1  → O   (degree 1–2)
    label 2  → N   (degree 1–3; degree-4 atoms handled as N⁺)
    label 3  → C   (degree 1–4)
    label 4  → S   (degree 1–4)
    label 5  → Cl  (mostly degree 1)
    label 6  → Br  (degree 1; degree-3/4 cases fall back to None)
    label 7  → F   (all degree 1)
    label 8  → I
    label 9  → P
    label 38 → Se
    labels 10–37 → trace elements (Na, K, Li, Ca, Sn, Mg, Si, As, B, …)

  Because bond orders are inferred from valence rules alone, the resulting
  SMILES represent the *topology* of each molecule exactly but may use a
  different Kekulé form or bond-order assignment than the original compound.
  Molecules whose sanitization cannot be fixed are silently skipped.

Output: data/processed/NCI109/nci109_smiles.smi
  Plain-text, space-separated, two columns:
      smiles zinc_id
  where zinc_id is the 1-based graph index from the TU dataset.
  Encoding: UTF-8, Unix line endings.

How to run (from the project root):
    python scripts/download/parse_nci109_to_smiles.py

Or with explicit paths:
    python scripts/download/parse_nci109_to_smiles.py \\
        --data-dir data/raw/NCI109/NCI109 \\
        --out data/processed/NCI109/nci109_smiles.smi
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT_DIR))

from rdkit import Chem
from rdkit.Chem import RWMol, SanitizeMol


ATOM_MAP: dict[int, str] = {
    1:  "O",
    2:  "N",
    3:  "C",
    4:  "S",
    5:  "Cl",
    6:  "Br",
    7:  "F",
    8:  "I",
    9:  "P",
    10: "Na",
    11: "K",
    12: "Li",
    13: "Ca",
    14: "Sn",
    15: "Mg",
    16: "Si",
    17: "As",
    18: "B",
    19: "Ge",
    20: "Pb",
    21: "Fe",
    22: "Zn",
    23: "Cu",
    24: "Mo",
    25: "Co",
    26: "Ag",
    27: "Al",
    28: "Ni",
    29: "Sb",
    30: "Cr",
    31: "Mn",
    32: "Hg",
    33: "Bi",
    34: "Au",
    35: "Ru",
    36: "Tl",
    37: "Te",
    38: "Se",
}

_NORMAL_VALENCE: dict[str, int] = {
    "C": 4, "N": 3, "O": 2, "S": 6, "F": 1, "Cl": 7, "Br": 7, "I": 7,
    "P": 5, "Na": 1, "K": 1, "Li": 1, "Ca": 2, "Sn": 4, "Mg": 2,
    "Si": 4, "As": 5, "B": 3, "Ge": 4, "Pb": 4, "Fe": 3, "Zn": 2,
    "Cu": 2, "Mo": 6, "Co": 3, "Ag": 1, "Al": 3, "Ni": 3, "Sb": 5,
    "Cr": 6, "Mn": 7, "Hg": 2, "Bi": 5, "Au": 3, "Ru": 8, "Tl": 3,
    "Te": 6, "Se": 6,
}



def _read_column(path: Path) -> list[int]:
    return [int(line.strip()) for line in path.read_text().splitlines() if line.strip()]


def _read_edges(path: Path) -> list[tuple[int, int]]:
    edges = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        u, v = [int(x) for x in line.split(",")]
        edges.append((u, v))
    return edges



def parse_nci109(data_dir: Path) -> list[dict]:
    graph_indicator = _read_column(data_dir / "NCI109_graph_indicator.txt")
    node_labels     = _read_column(data_dir / "NCI109_node_labels.txt")
    graph_labels    = _read_column(data_dir / "NCI109_graph_labels.txt")
    edges           = _read_edges(data_dir / "NCI109_A.txt")

    nodes_per_graph: dict[int, dict[int, str]] = defaultdict(dict)
    for node_id, (graph_id, label_idx) in enumerate(
        zip(graph_indicator, node_labels), start=1
    ):
        nodes_per_graph[graph_id][node_id] = ATOM_MAP[label_idx]

    edges_per_graph: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for u, v in edges:
        if u < v:
            graph_id = graph_indicator[u - 1]
            edges_per_graph[graph_id].append((u, v))

    graphs = []
    for graph_id, label in enumerate(graph_labels, start=1):
        graphs.append(
            {
                "graph_id": graph_id,
                "label":    label,
                "nodes":    nodes_per_graph[graph_id],
                "edges":    edges_per_graph[graph_id],
            }
        )
    return graphs



def _assign_formal_charges(mol: RWMol) -> None:
    """
    Heuristically fix formal charges for atoms that exceed their normal valence
    when all bonds are treated as single.

    Covers the common over-valent cases in NCI-family datasets:
      Cl/Br/I even degree  → +1  (neutral halogens only permit odd valences
                                   [1,3,5,7]; +1 shifts the list to [0,2,4,6])
      N  degree 4          → N⁺  (quaternary ammonium)
      O  degree 3          → O⁺  (oxonium)
      S  degree 3–4        → S⁺  (sulfonium)
      Se degree 3–4        → Se⁺ (selenonium)
    """
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if atom.GetFormalCharge() != 0:
            continue
        degree = atom.GetDegree()

        if sym in ("Cl", "Br", "I") and degree > 0 and degree % 2 == 0:
            atom.SetFormalCharge(1)
            continue

        if sym == "F" and degree > 1:
            atom.SetFormalCharge(degree - 1)
            continue

        normal = _NORMAL_VALENCE.get(sym)
        if normal is None or degree <= normal:
            continue
        if sym == "N" and degree == 4:
            atom.SetFormalCharge(1)
        elif sym == "O" and degree == 3:
            atom.SetFormalCharge(1)
        elif sym == "S" and degree in (3, 4):
            atom.SetFormalCharge(1)
        elif sym == "Se" and degree in (3, 4):
            atom.SetFormalCharge(1)



def graph_to_smiles(graph: dict) -> str | None:
    """
    Build an RDKit RWMol from node/edge data using single bonds throughout
    (NCI109 has no edge labels), apply formal-charge correction, sanitize,
    and return a canonical SMILES string.

    Returns None if the graph has no edges, or if both full and relaxed
    sanitization fail, or if the generated SMILES fails a round-trip check.
    """
    nodes: dict[int, str] = graph["nodes"]
    edges: list[tuple[int, int]] = graph["edges"]

    if not edges:
        return None

    sorted_node_ids = sorted(nodes.keys())
    node_to_idx: dict[int, int] = {nid: i for i, nid in enumerate(sorted_node_ids)}

    mol = RWMol()
    for nid in sorted_node_ids:
        mol.AddAtom(Chem.Atom(nodes[nid]))
    for u, v in edges:
        mol.AddBond(node_to_idx[u], node_to_idx[v], Chem.rdchem.BondType.SINGLE)

    _assign_formal_charges(mol)

    try:
        SanitizeMol(mol)
        smiles = Chem.MolToSmiles(mol)
    except Exception:
        try:
            SanitizeMol(
                mol,
                Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES,
            )
            smiles = Chem.MolToSmiles(mol)
        except Exception as exc:
            print(f"  [warn] graph_id={graph['graph_id']} sanitization failed: {exc}")
            return None

    if Chem.MolFromSmiles(smiles) is None:
        print(
            f"  [warn] graph_id={graph['graph_id']} SMILES failed round-trip "
            f"(over-valent atom): {smiles}"
        )
        return None

    return smiles



def main() -> None:
    parser = argparse.ArgumentParser(description="Convert NCI109 dataset to SMILES (.smi).")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=ROOT_DIR / "data" / "raw" / "NCI109" / "NCI109",
        help="Directory containing NCI109_*.txt files",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT_DIR / "data" / "processed" / "NCI109" / "nci109_smiles.smi",
        help="Output .smi path (space-separated: smiles zinc_id)",
    )
    args = parser.parse_args()

    data_dir: Path = args.data_dir
    out_path: Path = args.out

    print(f"Parsing NCI109 from: {data_dir}")
    graphs = parse_nci109(data_dir)
    print(f"Found {len(graphs)} graphs")

    rows: list[tuple[str, int]] = []
    failed = 0
    for i, g in enumerate(graphs, start=1):
        smiles = graph_to_smiles(g)
        if smiles is None:
            failed += 1
        else:
            rows.append((smiles, g["graph_id"]))
        if i % 500 == 0:
            print(f"  Processed {i}/{len(graphs)} …")

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
