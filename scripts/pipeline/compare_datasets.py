"""
Multi-dataset comparison pipeline for WL-identifiability analysis.

Runs the full WL-identifiability pipeline on each specified .smi dataset and
produces:
  - results/<dataset_name>/
        <dataset_name>.db   — per-molecule SQLite results
        summary.txt         — human-readable per-dataset summary
  - results/summary.csv     — one aggregated row per dataset

Usage (from project root):
    python scripts/compare_datasets.py \\
        --data data/processed/MUTAG/mutag_smiles.smi \\
               data/processed/NCI1/nci1_smiles.smi \\
               data/processed/NCI109/nci109_smiles.smi \\
               data/raw/ZINC/zinc250k.smi

Optional flags:
    --out        results/          # root output directory
    --sample     300               # molecules used to estimate WL steps
    --cap        50                # max WL iterations cap
    --workers    1                 # parallel worker processes
    --wl-steps   N                 # skip estimation, use this fixed value
"""
from __future__ import annotations

import argparse
import csv
import sqlite3
import sys
import time
import traceback
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

import polars as pl

from wl_identifiability.experiments import (
    estimate_fixed_wl_steps_from_dataframe,
    run_molecule_analysis_pipeline,
)
from scripts.database.db import open_db, insert_row as db_insert_row



def load_smi(path: Path) -> pl.DataFrame:
    """
    Read a space-separated .smi file with a header row (smiles zinc_id).

    Returns a Polars DataFrame with columns [smiles, zinc_id] (both str).
    Rows where smiles is the literal word 'smiles' or 'smi' are treated as
    header lines and skipped.
    """
    rows: list[dict] = []
    with path.open(encoding="utf-8") as fh:
        for i, raw in enumerate(fh):
            line = raw.strip()
            if not line:
                continue
            parts = line.split()
            smiles = parts[0]
            if smiles.lower() in ("smiles", "smi"):
                continue
            zinc_id = parts[1] if len(parts) > 1 else f"mol_{i}"
            rows.append({"smiles": smiles, "zinc_id": str(zinc_id)})

    return pl.DataFrame(rows)



def _safe_mean(series: pl.Series) -> float:
    """Return the mean of a numeric series, ignoring nulls. Returns 0.0 if empty."""
    dropped = series.drop_nulls()
    if len(dropped) == 0:
        return 0.0
    return float(dropped.mean())


def _safe_sum(series: pl.Series) -> int:
    """Return the sum of a boolean/int series cast to Int64, ignoring nulls."""
    if series.dtype in (pl.Boolean,):
        series = series.cast(pl.Int64)
    dropped = series.drop_nulls()
    if len(dropped) == 0:
        return 0
    return int(dropped.sum())


def aggregate_results(result_df: pl.DataFrame, n_failed: int) -> dict:
    """
    Compute dataset-level summary metrics from the pipeline result DataFrame.

    The result_df contains one row per molecule with columns produced by
    analyze_single_molecule (see wl_identifiability/experiments.py).
    """
    n_total = len(result_df)
    n_ok = max(n_total - n_failed, 0)

    ok_df = result_df.filter(pl.col("ok") == True) if "ok" in result_df.columns else result_df

    agg: dict = {
        "n_total": n_total,
        "n_ok": n_ok,
        "n_failed": n_failed,
        "avg_nodes": _safe_mean(ok_df["n_nodes"]) if "n_nodes" in ok_df.columns else 0.0,
        "avg_edges": _safe_mean(ok_df["n_edges"]) if "n_edges" in ok_df.columns else 0.0,
        "avg_wl_iters_atom": _safe_mean(ok_df["wl_iters_atom"]) if "wl_iters_atom" in ok_df.columns else 0.0,
        "avg_n_colors_atom": _safe_mean(ok_df["n_colors_atom"]) if "n_colors_atom" in ok_df.columns else 0.0,
        "avg_n_colors_top": _safe_mean(ok_df["n_colors_top"]) if "n_colors_top" in ok_df.columns else 0.0,
        "avg_color_ratio": _safe_mean(ok_df["color_ratio_atom_to_top"]) if "color_ratio_atom_to_top" in ok_df.columns else 0.0,
        "avg_flipped_edges": _safe_mean(ok_df["n_flipped_edges"]) if "n_flipped_edges" in ok_df.columns else 0.0,
        "n_top_bouquet_forest": _safe_sum(ok_df["top_bouquet_forest_verdict"]) if "top_bouquet_forest_verdict" in ok_df.columns else 0,
        "n_atom_bouquet_forest": _safe_sum(ok_df["atom_bouquet_forest_verdict"]) if "atom_bouquet_forest_verdict" in ok_df.columns else 0,
    }

    agg["has_top_bouquet_forest"]  = agg["n_top_bouquet_forest"] > 0
    agg["has_atom_bouquet_forest"] = agg["n_atom_bouquet_forest"] > 0
    return agg



def build_summary_text(
    dataset_name: str,
    smi_path: Path,
    wl_steps: int,
    k_p95: int,
    sample_size: int,
    total_runtime_s: float,
    agg: dict,
) -> str:
    """Return a human-readable summary string for one dataset run."""
    n_total = agg["n_total"]
    n_ok = agg["n_ok"]
    ok_pct = (100.0 * n_ok / n_total) if n_total > 0 else 0.0
    avg_rt = (total_runtime_s / n_ok) if n_ok > 0 else 0.0

    lines = [
        "=== WL-identifiability pipeline summary ===",
        f"Dataset           : {dataset_name}",
        f"Source file       : {smi_path}",
        f"Sample size (est) : {sample_size}",
        f"WL steps used     : {wl_steps}",
        f"K_p95             : {k_p95}",
        "",
        f"Total molecules   : {n_total:>8,}",
        f"Parsed OK         : {n_ok:>8,}  ({ok_pct:.1f}%)",
        f"Failed            : {agg['n_failed']:>8,}",
        "",
        f"Avg nodes         : {agg['avg_nodes']:>11.3f}",
        f"Avg edges         : {agg['avg_edges']:>11.3f}",
        f"Avg WL iters      : {agg['avg_wl_iters_atom']:>11.3f}",
        f"Avg colors (atom) : {agg['avg_n_colors_atom']:>11.3f}",
        f"Avg colors (top)  : {agg['avg_n_colors_top']:>11.3f}",
        f"Avg color ratio   : {agg['avg_color_ratio']:>11.3f}",
        f"Avg flipped edges : {agg['avg_flipped_edges']:>11.3f}",
        "",
        f"BF verdict (top)  : {agg['n_top_bouquet_forest']:>8,}",
        f"BF verdict (atom) : {agg['n_atom_bouquet_forest']:>8,}",
        "",
        f"Total runtime     : {total_runtime_s:>8.1f} s",
        f"Avg time/molecule : {avg_rt:>11.4f} s",
    ]
    return "\n".join(lines)



def finalize_db(db_path: Path) -> None:
    """Checkpoint WAL and switch back to DELETE journal mode so no -wal/-shm files remain."""
    if not db_path.exists():
        return
    with sqlite3.connect(str(db_path)) as conn:
        conn.execute("PRAGMA wal_checkpoint(FULL);")
        conn.execute("PRAGMA journal_mode=DELETE;")
        conn.commit()
    for suffix in ("-wal", "-shm"):
        extra = db_path.with_name(db_path.name + suffix)
        try:
            if extra.exists():
                extra.unlink()
        except PermissionError:
            print(f"  [WARN] Could not remove {extra.name} — file still open.")



def run_one_dataset(
    smi_path: Path,
    out_root: Path,
    sample_size: int,
    cap: int,
    workers: int,
    fixed_wl_steps: int | None,
) -> dict | None:
    """
    Run the full pipeline on one .smi file.

    Returns a dict of aggregated metrics (suitable for summary.csv) or None if
    the dataset could not be loaded or the pipeline raised an unhandled error.
    """
    name = smi_path.stem
    print(f"\n{'='*60}")
    print(f"  Dataset : {name}")
    print(f"  File    : {smi_path}")
    print(f"{'='*60}")

    try:
        df = load_smi(smi_path)
    except Exception:
        print(f"  [ERROR] Failed to load {smi_path}:")
        traceback.print_exc()
        return None

    if len(df) == 0:
        print(f"  [ERROR] No molecules loaded from {smi_path}. Skipping.")
        return None

    print(f"  Loaded {len(df):,} molecules.")

    if fixed_wl_steps is not None:
        wl_steps = fixed_wl_steps
        k_p95 = fixed_wl_steps
        print(f"  WL steps: {wl_steps} (user-specified)")
    else:
        print(f"  Estimating WL steps (sample={min(sample_size, len(df))}, cap={cap}) …")
        k_stats = estimate_fixed_wl_steps_from_dataframe(df, sample_size=sample_size, cap=cap)
        wl_steps = int(k_stats["K_max"])
        k_p95 = int(k_stats["K_p95"])
        print(f"  WL steps: {wl_steps}  (K_p95={k_p95})")

    run_dir = out_root / name
    run_dir.mkdir(parents=True, exist_ok=True)
    db_path = run_dir / f"{name}.db"

    print(f"  Running pipeline on {len(df):,} molecules (workers={workers}) …")
    t0 = time.perf_counter()
    try:
        result_df, n_failed = run_molecule_analysis_pipeline(
            df, fixed_wl_steps=wl_steps, n_workers=workers
        )
    except Exception:
        print(f"  [ERROR] Pipeline crashed for dataset '{name}':")
        traceback.print_exc()
        return None
    elapsed = time.perf_counter() - t0
    print(f"  Done in {elapsed:.1f}s  ({n_failed} failed)")

    conn = open_db(db_path)
    for row in result_df.iter_rows(named=True):
        try:
            db_insert_row(conn, row)
        except Exception as exc:
            print(f"  [WARN] DB insert failed for {row.get('zinc_id')}: {exc}")
    conn.close()
    finalize_db(db_path)

    agg = aggregate_results(result_df, n_failed)
    avg_rt = elapsed / agg["n_ok"] if agg["n_ok"] > 0 else 0.0

    summary_text = build_summary_text(
        dataset_name=name,
        smi_path=smi_path,
        wl_steps=wl_steps,
        k_p95=k_p95,
        sample_size=min(sample_size, len(df)),
        total_runtime_s=elapsed,
        agg=agg,
    )
    print("\n" + summary_text)

    summary_path = run_dir / "summary.txt"
    summary_path.write_text(summary_text, encoding="utf-8")
    print(f"\n  Summary : {summary_path}")
    print(f"  SQLite  : {db_path}")

    return {
        "dataset": name,
        "smi_file": str(smi_path),
        "n_total": agg["n_total"],
        "n_ok": agg["n_ok"],
        "n_failed": agg["n_failed"],
        "avg_nodes": round(agg["avg_nodes"], 4),
        "avg_edges": round(agg["avg_edges"], 4),
        "avg_wl_iters_atom": round(agg["avg_wl_iters_atom"], 4),
        "avg_n_colors_atom": round(agg["avg_n_colors_atom"], 4),
        "avg_n_colors_top": round(agg["avg_n_colors_top"], 4),
        "avg_color_ratio": round(agg["avg_color_ratio"], 4),
        "avg_flipped_edges": round(agg["avg_flipped_edges"], 4),
        "n_top_bouquet_forest": agg["n_top_bouquet_forest"],
        "n_atom_bouquet_forest": agg["n_atom_bouquet_forest"],
        "has_top_bouquet_forest": agg["has_top_bouquet_forest"],
        "has_atom_bouquet_forest": agg["has_atom_bouquet_forest"],
        "wl_steps_used": wl_steps,
        "k_p95": k_p95,
        "total_runtime_s": round(elapsed, 3),
        "avg_runtime_per_mol_s": round(avg_rt, 6),
    }



def resolve_smi_paths(inputs: list[str], recursive: bool = False) -> list[Path]:
    """
    Accept a mixed list of .smi file paths and/or directories.

    For each entry:
      - If it is an existing .smi file, include it directly.
      - If it is a directory, glob for all *.smi files inside it
        (recursively when recursive=True, otherwise top-level only).

    The final list is de-duplicated and sorted by name for deterministic runs.
    """
    collected: list[Path] = []

    for raw in inputs:
        p = Path(raw)
        if p.is_dir():
            pattern = "**/*.smi" if recursive else "*.smi"
            found = sorted(p.glob(pattern), key=lambda x: x.name)
            if not found:
                print(f"[WARN] No .smi files found in directory: {p}")
            collected.extend(found)
        elif p.suffix.lower() == ".smi":
            collected.append(p)
        else:
            collected.append(p)

    seen: set[Path] = set()
    unique: list[Path] = []
    for p in collected:
        resolved = p.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique.append(p)

    return unique



CSV_FIELDS = [
    "dataset",
    "smi_file",
    "n_total",
    "n_ok",
    "n_failed",
    "avg_nodes",
    "avg_edges",
    "avg_wl_iters_atom",
    "avg_n_colors_atom",
    "avg_n_colors_top",
    "avg_color_ratio",
    "avg_flipped_edges",
    "n_top_bouquet_forest",
    "n_atom_bouquet_forest",
    "has_top_bouquet_forest",
    "has_atom_bouquet_forest",
    "wl_steps_used",
    "k_p95",
    "total_runtime_s",
    "avg_runtime_per_mol_s",
]


def write_summary_csv(rows: list[dict], csv_path: Path) -> None:
    """Write aggregated per-dataset rows to a UTF-8 CSV file."""
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", encoding="utf-8", newline="\n") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)



def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run the WL-identifiability pipeline on multiple datasets and compare results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument(
        "--data",
        nargs="+",
        required=True,
        metavar="PATH",
        help=(
            "One or more .smi files and/or directories. "
            "When a directory is given, all *.smi files inside it are used."
        ),
    )
    p.add_argument(
        "--recursive",
        action="store_true",
        help="When --data includes a directory, search for .smi files recursively.",
    )
    p.add_argument(
        "--out",
        default="results",
        metavar="DIR",
        help="Root output directory (default: results/).",
    )
    p.add_argument(
        "--sample",
        type=int,
        default=300,
        help="Molecules sampled for WL-step estimation (default: 300).",
    )
    p.add_argument(
        "--cap",
        type=int,
        default=50,
        help="Maximum WL iterations cap (default: 50).",
    )
    p.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Parallel worker processes (default: 1).",
    )
    p.add_argument(
        "--wl-steps",
        type=int,
        default=None,
        metavar="N",
        help="Fixed WL step count for all datasets (skips per-dataset estimation).",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)

    smi_paths = resolve_smi_paths(args.data, recursive=getattr(args, "recursive", False))
    if not smi_paths:
        print("[ERROR] No .smi files found for the given --data inputs.")
        return 1

    print(f"Found {len(smi_paths)} .smi file(s) to process:")
    for p in smi_paths:
        print(f"  {p}")

    summary_rows: list[dict] = []
    failed_datasets: list[str] = []

    for smi_path in smi_paths:
        if not smi_path.exists():
            print(f"\n[ERROR] File not found: {smi_path}. Skipping.")
            failed_datasets.append(str(smi_path))
            continue

        row = run_one_dataset(
            smi_path=smi_path,
            out_root=out_root,
            sample_size=args.sample,
            cap=args.cap,
            workers=args.workers,
            fixed_wl_steps=args.wl_steps,
        )

        if row is None:
            failed_datasets.append(smi_path.stem)
        else:
            summary_rows.append(row)

    if summary_rows:
        csv_path = out_root / "summary.csv"
        write_summary_csv(summary_rows, csv_path)
        print(f"\n{'='*60}")
        print(f"Global summary CSV : {csv_path}")
        print(f"Datasets processed : {len(summary_rows)}")
        if failed_datasets:
            print(f"Datasets failed    : {', '.join(failed_datasets)}")
    else:
        print("\n[ERROR] No datasets were processed successfully.")
        return 1

    return 0 if not failed_datasets else 1


if __name__ == "__main__":
    sys.exit(main())
