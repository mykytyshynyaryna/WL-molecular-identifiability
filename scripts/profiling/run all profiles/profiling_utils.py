"""
Shared profiling utilities.

Provides:
  - timed_run()     : measure wall-clock time over N repetitions
  - cprofile_run()  : run cProfile over a function, return Stats object
  - print_stats()   : print top-N slowest functions from a Stats object
  - save_stats()    : save .prof binary and .txt report to disk
  - build_table()   : render a list-of-dicts as an ASCII table
"""
from __future__ import annotations

import cProfile
import io
import pstats
import time
from pathlib import Path


def timed_run(fn, cases: list, reps: int = 1) -> float:
    """
    Run fn(G, labels) for every case, repeated reps times.
    Returns total elapsed wall-clock seconds.
    """
    t0 = time.perf_counter()
    for _ in range(reps):
        for _, G, labels in cases:
            fn(G, labels)
    return time.perf_counter() - t0



def cprofile_run(fn, cases: list) -> pstats.Stats:
    """
    Run cProfile over one pass of fn(G, labels) for every case.
    Returns a pstats.Stats object.
    """
    pr = cProfile.Profile()
    pr.enable()
    for _, G, labels in cases:
        fn(G, labels)
    pr.disable()

    stream = io.StringIO()
    stats = pstats.Stats(pr, stream=stream)
    stats.sort_stats("cumulative")
    return stats


def print_stats(stats: pstats.Stats, top_n: int = 20) -> None:
    """Print the top_n slowest functions from a Stats object."""
    stream = io.StringIO()
    stats.stream = stream
    stats.print_stats(top_n)
    print(stream.getvalue())


def save_stats(
    stats: pstats.Stats,
    out_dir: Path,
    stem: str,
    top_n: int = 40,
) -> None:
    """
    Save profiling results to out_dir/<stem>.prof and out_dir/<stem>.txt.
    Creates out_dir if it does not exist.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    prof_path = out_dir / f"{stem}.prof"
    stats.dump_stats(str(prof_path))

    txt_path = out_dir / f"{stem}.txt"
    stream = io.StringIO()
    stats.stream = stream
    stats.print_stats(top_n)
    txt_path.write_text(stream.getvalue(), encoding="utf-8")

    print(f"  Saved .prof → {prof_path}")
    print(f"  Saved .txt  → {txt_path}")



def build_table(rows: list[dict]) -> str:
    """Render a list of dicts as a fixed-width ASCII table."""
    if not rows:
        return ""
    fields = list(rows[0].keys())
    col_w = {
        f: max(len(f), max(len(str(r[f])) for r in rows))
        for f in fields
    }
    sep    = "+" + "+".join("-" * (col_w[f] + 2) for f in fields) + "+"
    header = "|" + "|".join(f" {f:<{col_w[f]}} " for f in fields) + "|"
    lines  = [sep, header, sep]
    for r in rows:
        lines.append("|" + "|".join(f" {str(r[f]):<{col_w[f]}} " for f in fields) + "|")
    lines.append(sep)
    return "\n".join(lines)
