"""
Profile: flip graph construction — NEW (optimized NumPy submatrix) implementation.

Function profiled : build_flip_graph_new (defined in this file)
Implementation    : same submatrix-flip scheme as the old version, but with:
                    - int8 adjacency matrix (8x smaller than float64, better cache)
                    - matrix built directly from the edge list (no nx.to_numpy_array)
                    - color-class node IDs pre-converted to numpy index arrays once
                    - integer threshold comparison (no float division)
                    - correct node relabeling for non-0-indexed graphs
Compared against  : profile_flip_graph_old.py (float64, nx.to_numpy_array, list indices)

What is measured:
  - cProfile breakdown of a single pass over the dataset
  - Total wall-clock time over --reps repetitions
  - Top slowest functions by cumulative time

Usage (from project root):
    python profiling/profile_flip_graph_new.py
    python profiling/profile_flip_graph_new.py --data data/AAAC.smi --max-mols 500
    python profiling/profile_flip_graph_new.py --reps 10 --no-save

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
import sys
import numpy as np
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).parent.parent))

from profiling.dataset_loader import load_cases
from profiling.profiling_utils import (
    build_table,
    cprofile_run,
    print_stats,
    save_stats,
    timed_run,
)

LABEL = "flip_graph_new"



def _group_nodes_by_color(labels: dict) -> dict:
    color2nodes: dict = {}
    for v, c in labels.items():
        color2nodes.setdefault(c, []).append(v)
    return color2nodes


def build_flip_graph_new(G: nx.Graph, labels: dict) -> tuple[nx.Graph, dict]:
    nodes = list(G.nodes())
    n = len(nodes)
    node_to_idx = {v: i for i, v in enumerate(nodes)}

    A = np.zeros((n, n), dtype=np.int8)
    for u, v in G.edges():
        iu, iv = node_to_idx[u], node_to_idx[v]
        A[iu, iv] = A[iv, iu] = 1

    color2nodes = _group_nodes_by_color(labels)
    colors = list(color2nodes.keys())

    color_idx = [
        np.array([node_to_idx[v] for v in color2nodes[c]], dtype=np.intp)
        for c in colors
    ]

    info = {"within_copy": 0, "within_flip": 0, "between_copy": 0, "between_flip": 0}

    for idx in color_idx:
        ni = len(idx)
        if ni < 2:
            continue
        ix = np.ix_(idx, idx)
        sub = A[ix]
        if sub.sum() > (ni * (ni - 1)) >> 1:
            flip = 1 - sub
            np.fill_diagonal(flip, 0)
            A[ix] = flip
            info["within_flip"] += 1
        else:
            info["within_copy"] += 1

    for ci in range(len(colors)):
        for cj in range(ci + 1, len(colors)):
            idx_i = color_idx[ci]
            idx_j = color_idx[cj]
            ni, nj = len(idx_i), len(idx_j)
            if ni > nj:
                idx_i, idx_j = idx_j, idx_i
                ni, nj = nj, ni
            ix_ij = np.ix_(idx_i, idx_j)
            sub = A[ix_ij]
            if sub.sum() * 2 > ni * nj:
                flip = 1 - sub
                A[ix_ij] = flip
                A[np.ix_(idx_j, idx_i)] = flip.T
                info["between_flip"] += 1
            else:
                info["between_copy"] += 1

    F = nx.from_numpy_array(A)
    F = nx.relabel_nodes(F, {i: v for i, v in enumerate(nodes)})
    return F, info


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
    stats = cprofile_run(build_flip_graph_new, cases)
    print_stats(stats, top_n=args.top_n)

    print(f"=== Timing ({args.reps} reps × {len(cases)} molecules) ===")
    total_s = timed_run(build_flip_graph_new, cases, reps=args.reps)
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
