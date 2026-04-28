"""
Profile: flip graph construction — OLD (NumPy submatrix) implementation.

Function profiled : build_flip_graph (defined in this file)
Implementation    : converts the graph to a dense NumPy adjacency matrix once,
                    then edits submatrix blocks in-place using np.ix_ indexing.
                    No Python edge-iteration loops.
Compared against  : profile_flip_graph_new.py (pure edge-list, no dense matrix)

What is measured:
  - cProfile breakdown of a single pass over the dataset
  - Total wall-clock time over --reps repetitions
  - Top slowest functions by cumulative time

Usage (from project root):
    python profiling/profile_flip_graph_old.py
    python profiling/profile_flip_graph_old.py --data data/AAAC.smi --max-mols 500
    python profiling/profile_flip_graph_old.py --reps 10 --no-save

Arguments:
    --data      Path to .smi dataset file          (default: data/AAAA.smi)
    --max-mols  Maximum molecules to load          (default: all)
    --reps      Timing repetitions                 (default: 5)
    --out       Output directory for saved reports (default: profiling/results)
    --no-save   Skip saving .prof/.txt files
    --top-n     Number of slowest functions shown  (default: 20)
"""
from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from profiling.dataset_loader import load_cases
from profiling.profiling_utils import (
    build_table,
    cprofile_run,
    print_stats,
    save_stats,
    timed_run,
)

LABEL = "flip_graph_old"



def _group_nodes_by_color(labels: dict) -> dict:
    color2nodes: dict = {}
    for v, c in labels.items():
        color2nodes.setdefault(c, []).append(v)
    return color2nodes


def build_flip_graph_old(G: nx.Graph, labels: dict) -> nx.Graph:
    A = nx.to_numpy_array(G)
    color_nodes = _group_nodes_by_color(labels)

    for i in color_nodes.values():
        subN = len(i)
        if subN < 2:
            continue
        subA = A[np.ix_(i, i)]
        if subA.sum() > (subN * (subN - 1)) / 2:
            flip_matrix = 1 - subA
            np.fill_diagonal(flip_matrix, 0)
            A[np.ix_(i, i)] = flip_matrix

    for i, j in itertools.combinations(color_nodes.values(), 2):
        subB = A[np.ix_(i, j)]
        if subB.sum() > (len(i) * len(j)) / 2:
            flip_matrix = 1 - subB
            A[np.ix_(i, j)] = flip_matrix
            A[np.ix_(j, i)] = flip_matrix.T

    return nx.from_numpy_array(A)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data",     default="data/AAAA.smi")
    parser.add_argument("--max-mols", type=int, default=None)
    parser.add_argument("--reps",     type=int, default=5)
    parser.add_argument("--out",      default="profiling/results")
    parser.add_argument("--no-save",  action="store_true")
    parser.add_argument("--top-n",    type=int, default=20)
    args = parser.parse_args()

    data_path = Path(args.data)
    if not data_path.exists():
        print(f"ERROR: file not found: {data_path}")
        sys.exit(1)

    print(f"Loading molecules from {data_path} ...")
    cases = load_cases(data_path, max_mols=args.max_mols)
    print(f"Loaded {len(cases)} molecules.\n")

    if not cases:
        print("No valid molecules found.")
        sys.exit(1)

    print("=== cProfile (1 pass) ===")
    stats = cprofile_run(build_flip_graph_old, cases)
    print_stats(stats, top_n=args.top_n)

    print(f"=== Timing ({args.reps} reps × {len(cases)} molecules) ===")
    total_s = timed_run(build_flip_graph_old, cases, reps=args.reps)
    n_calls = len(cases) * args.reps
    rows = [
        {"metric": "total_s",  "value": f"{total_s:.4f}"},
        {"metric": "ms/call",  "value": f"{1000 * total_s / n_calls:.4f}"},
        {"metric": "reps",     "value": str(args.reps)},
        {"metric": "n_mols",   "value": str(len(cases))},
    ]
    print(build_table(rows))

    if not args.no_save:
        save_stats(stats, Path(args.out), LABEL, top_n=40)


if __name__ == "__main__":
    main()
