"""
Profile: bouquet detection — BASELINE implementation.

Function profiled : wl_identifiability.bouquet._check_bouquet_component_baseline
Implementation    : uses nx.cycle_basis, graph copies, and nx.is_isomorphic for
                    pairwise petal comparison.
Compared against  : profile_bouquet_optimized.py (leaf-stripping + AHU signatures)

What is measured:
  - cProfile breakdown of a single pass over the dataset
  - Total wall-clock time over --reps repetitions
  - Top slowest functions by cumulative time

Usage (from project root):
    python profiling/profile_bouquet_baseline.py
    python profiling/profile_bouquet_baseline.py --data data/AAAC.smi --max-mols 500
    python profiling/profile_bouquet_baseline.py --reps 10 --no-save

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
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from wl_identifiability.bouquet import _check_bouquet_component_baseline
from profiling.dataset_loader import load_cases
from profiling.profiling_utils import (
    build_table,
    cprofile_run,
    print_stats,
    save_stats,
    timed_run,
)

LABEL = "bouquet_baseline"


def _run(G, labels):
    return _check_bouquet_component_baseline(G, labels=labels)


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
    stats = cprofile_run(_run, cases)
    print_stats(stats, top_n=args.top_n)

    print(f"=== Timing ({args.reps} reps × {len(cases)} molecules) ===")
    total_s = timed_run(_run, cases, reps=args.reps)
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
