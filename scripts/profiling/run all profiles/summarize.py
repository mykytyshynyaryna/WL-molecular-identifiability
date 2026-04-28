"""
Parse all profiling result .txt files in --results-dir and print a
comparative summary table.

Extracts from each cProfile report:
  - total function calls
  - total runtime (seconds)
  - top-3 hottest functions

Computes speedup between paired scenarios:
  - bouquet_optimized  vs  bouquet_baseline
  - flip_graph_new     vs  flip_graph_old

Usage (from project root):
    python profiling/summarize.py
    python profiling/summarize.py --results-dir profiling/results
    python profiling/summarize.py --n-mols 3730 --no-save
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from profiling.profiling_utils import build_table

SCENARIO_ORDER = [
    "bouquet_baseline",
    "bouquet_optimized",
    "flip_graph_old",
    "flip_graph_new",
]

LABELS = {
    "bouquet_baseline":  "bouquet  baseline  (_check_bouquet_component_baseline)",
    "bouquet_optimized": "bouquet  optimized (_check_bouquet_component_optimized)",
    "flip_graph_old":    "flip graph  old    (NumPy submatrix, float64)",
    "flip_graph_new":    "flip graph  new    (NumPy submatrix, int8, pre-indexed)",
}

SPEEDUP_PAIRS = [
    ("bouquet_optimized",  "bouquet_baseline"),
    ("flip_graph_new",     "flip_graph_old"),
]



_RE_HEADER = re.compile(
    r"(\d+)\s+function calls.*in\s+([\d.]+)\s+seconds"
)
_RE_ROW = re.compile(
    r"^\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(.+)$"
)


def _basename(path_str: str) -> str:
    """Extract bare filename from a full path string inside a cProfile line."""
    p = path_str.strip()
    m = re.search(r"\(([^)]+)\)\s*$", p)
    return m.group(1) if m else p


def parse_profile_txt(path: Path) -> dict:
    """
    Parse a cProfile .txt report.

    Returns a dict with keys:
        total_calls, total_s, top_functions (list of dicts)
    """
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    total_calls = 0
    total_s = 0.0
    top_functions: list[dict] = []

    for line in lines:
        m = _RE_HEADER.search(line)
        if m:
            total_calls = int(m.group(1))
            total_s = float(m.group(2))
            continue

        m = _RE_ROW.match(line)
        if m and len(top_functions) < 3:
            ncalls   = m.group(1)
            tottime  = float(m.group(2))
            cumtime  = float(m.group(4))
            location = m.group(6)
            top_functions.append({
                "ncalls":   ncalls,
                "tottime":  tottime,
                "cumtime":  cumtime,
                "function": _basename(location),
            })

    return {
        "total_calls": total_calls,
        "total_s":     total_s,
        "top":         top_functions,
    }



def build_summary(results: dict[str, dict], n_mols: int) -> str:
    """
    Build the full summary text from parsed results.

    results : {scenario_name: parse_profile_txt(...)}
    n_mols  : number of molecules profiled (for ms/call calculation)
    """
    lines: list[str] = []

    lines.append("=" * 72)
    lines.append("  PROFILING SUMMARY")
    lines.append(f"  molecules: {n_mols:,}")
    lines.append("=" * 72)
    lines.append("")

    overview_rows = []
    for name in SCENARIO_ORDER:
        if name not in results:
            continue
        r = results[name]
        ms = (r["total_s"] / n_mols * 1000) if n_mols > 0 else 0.0
        overview_rows.append({
            "scenario":       name,
            "total_s":        f"{r['total_s']:.3f}",
            "ms/mol":         f"{ms:.4f}",
            "total_calls":    f"{r['total_calls']:,}",
        })

    lines.append("  Overview (1 pass, all molecules)")
    lines.append("")
    lines.append(build_table(overview_rows))
    lines.append("")

    speedup_rows = []
    for faster, slower in SPEEDUP_PAIRS:
        if faster not in results or slower not in results:
            continue
        t_fast = results[faster]["total_s"]
        t_slow = results[slower]["total_s"]
        c_fast = results[faster]["total_calls"]
        c_slow = results[slower]["total_calls"]
        speedup_rows.append({
            "comparison":      f"{faster}  vs  {slower}",
            "time_speedup":    f"{t_slow / t_fast:.2f}x  faster",
            "calls_reduction": f"{c_slow / c_fast:.2f}x  fewer calls",
        })

    if speedup_rows:
        lines.append("  Speedup comparison")
        lines.append("")
        lines.append(build_table(speedup_rows))
        lines.append("")

    lines.append("  Top-3 hottest functions per scenario")
    lines.append("")

    for name in SCENARIO_ORDER:
        if name not in results:
            continue
        r = results[name]
        label = LABELS.get(name, name)
        lines.append(f"  [{name}]  {label}")
        lines.append(f"  total: {r['total_s']:.3f}s   calls: {r['total_calls']:,}")

        top_rows = [
            {
                "rank":     str(i + 1),
                "ncalls":   fn["ncalls"],
                "tottime":  f"{fn['tottime']:.3f}",
                "cumtime":  f"{fn['cumtime']:.3f}",
                "function": fn["function"][:55],
            }
            for i, fn in enumerate(r["top"])
        ]
        if top_rows:
            lines.append(build_table(top_rows))
        lines.append("")

    lines.append("=" * 72)
    return "\n".join(lines)



def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--results-dir", default="profiling/results",
        help="Directory containing cProfile .txt reports",
    )
    parser.add_argument(
        "--n-mols", type=int, default=3730,
        help="Number of molecules in the profiled dataset (for ms/mol). "
             "Pass the total count when --data-dir was used.",
    )
    parser.add_argument(
        "--no-save", action="store_true",
        help="Print to stdout only, do not write summary.txt",
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    if not results_dir.exists():
        print(f"ERROR: results directory not found: {results_dir}")
        sys.exit(1)

    results: dict[str, dict] = {}
    for name in SCENARIO_ORDER:
        txt_path = results_dir / f"{name}.txt"
        if txt_path.exists():
            results[name] = parse_profile_txt(txt_path)
        else:
            print(f"  [WARN] {txt_path} not found — skipping")

    if not results:
        print("No profile result files found.")
        sys.exit(1)

    summary = build_summary(results, n_mols=args.n_mols)
    print(summary)

    if not args.no_save:
        out_path = results_dir / "summary.txt"
        out_path.write_text(summary, encoding="utf-8")
        print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
