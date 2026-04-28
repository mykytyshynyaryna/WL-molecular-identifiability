"""
Run all profiling scenarios sequentially and print a comparative summary.

Each scenario is run with the same dataset, molecule count, and repetitions.
Results are saved to --out (default: profiling/results).

Usage (from project root):
    python profiling/run_all_profiles.py
    python profiling/run_all_profiles.py --data data/AAAC.smi --max-mols 500
    python profiling/run_all_profiles.py --data-dir data/              # all .smi files
    python profiling/run_all_profiles.py --reps 10 --no-save

Arguments:
    --data      Path to a single .smi dataset file (default: data/AAAA.smi)
    --data-dir  Directory with .smi files — loads and combines all of them
    --max-mols  Maximum molecules to load total   (default: all)
    --reps      Timing repetitions per scenario   (default: 5)
    --out       Output directory for saved reports (default: profiling/results)
    --no-save   Skip saving .prof/.txt files

Note: --data-dir takes precedence over --data when both are provided.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import itertools

import networkx as nx
import numpy as np

from wl_identifiability.bouquet import (
    _check_bouquet_component_baseline,
    _check_bouquet_component_optimized,
)
from wl_identifiability.flip_graph import build_flip_graph_from_labels
from profiling.dataset_loader import load_cases, load_cases_from_dir
from profiling.profiling_utils import (
    build_table,
    cprofile_run,
    save_stats,
    timed_run,
)



def _group_nodes_by_color(labels: dict) -> dict:
    color2nodes: dict = {}
    for v, c in labels.items():
        color2nodes.setdefault(c, []).append(v)
    return color2nodes


def _build_flip_graph_old(G: nx.Graph, labels: dict) -> nx.Graph:
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



SCENARIOS = [
    ("bouquet_baseline",  lambda G, lbl: _check_bouquet_component_baseline(G, labels=lbl)),
    ("bouquet_optimized", lambda G, lbl: _check_bouquet_component_optimized(G, labels=lbl)),
    ("flip_graph_old",    _build_flip_graph_old),
    ("flip_graph_new",    lambda G, lbl: build_flip_graph_from_labels(G, lbl)),
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data",     default="data/AAAA.smi")
    parser.add_argument("--data-dir", default=None,
                        help="Directory with .smi files (overrides --data)")
    parser.add_argument("--max-mols", type=int, default=None)
    parser.add_argument("--reps",     type=int, default=5)
    parser.add_argument("--out",      default="profiling/results")
    parser.add_argument("--no-save",  action="store_true")
    args = parser.parse_args()

    if args.data_dir is not None:
        data_dir = Path(args.data_dir)
        if not data_dir.is_dir():
            print(f"ERROR: not a directory: {data_dir}")
            sys.exit(1)
        print(f"Loading molecules from all .smi files in {data_dir} ...")
        cases, n_files = load_cases_from_dir(data_dir, max_mols=args.max_mols)
        print(f"Loaded {len(cases)} molecules from {n_files} files.\n")
    else:
        data_path = Path(args.data)
        if not data_path.exists():
            print(f"ERROR: file not found: {data_path}")
            sys.exit(1)
        print(f"Loading molecules from {data_path} ...")
        cases = load_cases(data_path, max_mols=args.max_mols)
        n_files = 1
        print(f"Loaded {len(cases)} molecules.\n")

    if not cases:
        print("No valid molecules found.")
        sys.exit(1)

    summary_rows = []

    for name, fn in SCENARIOS:
        print(f"--- {name} ---")
        stats = cprofile_run(fn, cases)
        total_s = timed_run(fn, cases, reps=args.reps)
        n_calls = len(cases) * args.reps
        ms_per_call = 1000 * total_s / n_calls

        print(f"  total_s={total_s:.4f}   ms/call={ms_per_call:.4f}\n")

        if not args.no_save:
            save_stats(stats, Path(args.out), name, top_n=40)

        summary_rows.append({
            "scenario":  name,
            "total_s":   f"{total_s:.4f}",
            "ms/call":   f"{ms_per_call:.4f}",
            "reps":      str(args.reps),
            "n_mols":    str(len(cases)),
            "n_files":   str(n_files),
        })

    print("\n=== Summary ===")
    print(build_table(summary_rows))


if __name__ == "__main__":
    main()
