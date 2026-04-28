"""
Run the full WL-identifiability analysis pipeline on a .smi file.

Loads SMILES from a .smi file, runs the full pipeline (WL coloring,
flip graph, skeleton, bouquet forest), and saves results to a shared
SQLite database with a per-dataset summary row appended to a CSV file.

All molecule results go into one global DB (default: results/all_results.db),
tagged with the source dataset file name in the ``dataset_name`` column.
One summary row per dataset is appended to results/summary.csv
(existing rows for the same dataset are replaced on rerun).

Usage (from project root):
    python scripts/pipeline/run_pipeline.py --smi data/raw/ZINC/zinc250k.smi
    python scripts/pipeline/run_pipeline.py --smi data/raw/ZINC/zinc250k.smi --workers 4
    python scripts/pipeline/run_pipeline.py --smi data/raw/ZINC/zinc250k.smi --limit 1000
    python scripts/pipeline/run_pipeline.py --smi data/raw/ZINC/zinc250k.smi --out results/all_results.db
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT_DIR / "src"))  # src/ before root so wl_identifiability resolves to src/wl_identifiability/
sys.path.insert(1, str(ROOT_DIR))

import polars as pl

from wl_identifiability.experiments import (
    estimate_fixed_wl_steps_from_dataframe,
    run_molecule_analysis_pipeline,
)
from scripts.database.db import open_db, insert_row, write_dataset_summary_row


DEFAULT_OUT = ROOT_DIR / "results" / "all_results.db"
DEFAULT_SUMMARY_CSV = ROOT_DIR / "results" / "summary.csv"



def load_smi(path: Path, limit: int | None) -> pl.DataFrame:
    """Read a .smi file and return a DataFrame with columns [smiles, zinc_id]."""
    rows = []
    with path.open(encoding="utf-8") as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            smiles = parts[0]
            if smiles.lower() in ("smiles", "smi"):
                continue
            zinc_id = parts[1] if len(parts) > 1 else f"mol_{i}"
            rows.append({"smiles": smiles, "zinc_id": zinc_id})
            if limit is not None and len(rows) >= limit:
                break

    print(f"Loaded {len(rows):,} molecules from {path}")
    return pl.DataFrame(rows)



def _build_summary_row(
    result_df: pl.DataFrame,
    n_failed: int,
    dataset_name: str,
    sample_size: int,
    wl_cap: int,
    k_max: int,
    k_p95: int,
    wl_steps_used: int,
) -> dict:
    """Compute one summary row dict from the pipeline result DataFrame."""
    total = len(result_df)
    n_ok = max(total - n_failed, 0)
    ok_pct = round(100.0 * n_ok / total, 1) if total > 0 else 0.0

    def _col_sum(col: str) -> int:
        if col not in result_df.columns:
            return 0
        v = result_df[col].cast(pl.Int64, strict=False).sum()
        return int(v) if v is not None else 0

    def _col_mean(col: str) -> float:
        if col not in result_df.columns:
            return 0.0
        v = result_df[col].mean()
        return round(float(v), 3) if v is not None else 0.0

    def _col_max(col: str) -> int:
        if col not in result_df.columns:
            return 0
        v = result_df[col].max()
        return int(v) if v is not None else 0

    return {
        "dataset_name":      dataset_name,
        "sample_size":       sample_size,
        "wl_cap":            wl_cap,
        "k_max":             k_max,
        "k_p95":             k_p95,
        "wl_steps_used":     wl_steps_used,
        "total_molecules":   total,
        "parsed_ok":         n_ok,
        "parsed_ok_pct":     ok_pct,
        "parse_wl_errors":   n_failed,
        "bf_verdict_top":    _col_sum("top_bouquet_forest_verdict"),
        "bf_verdict_atom":   _col_sum("atom_bouquet_forest_verdict"),
        "avg_colors_top":    _col_mean("n_colors_top"),
        "avg_colors_atom":   _col_mean("n_colors_atom"),
        "avg_color_ratio":   _col_mean("color_ratio_atom_to_top"),
        "avg_flipped_edges": _col_mean("n_flipped_edges"),
        "max_flipped_edges": _col_max("n_flipped_edges"),
    }


def _build_summary_text(
    dataset_name: str,
    elapsed: float,
    n_total: int,
    n_ok: int,
    n_failed: int,
    out_path: Path,
    summary_csv: str,
    result_df: pl.DataFrame,
) -> str:
    """Build the human-readable per-dataset summary block (terminal output only)."""
    lines: list[str] = []
    lines.append(f"Done in {elapsed:.1f}s")
    lines.append(f"  Dataset  : {dataset_name}")
    lines.append(f"  Total    : {n_total:,}")
    lines.append(f"  OK       : {n_ok:,}")
    lines.append(f"  Failed   : {n_failed:,}")
    lines.append(f"  Results  : {out_path}")
    lines.append(f"  Summary  : {Path(summary_csv)}")

    for mode in ("top", "atom"):
        col = f"{mode}_bouquet_forest_verdict"
        if col in result_df.columns:
            verdicts = result_df[col].drop_nulls()
            n_true  = int((verdicts == 1).sum())
            n_false = int((verdicts == 0).sum())
            lines.append("")
            lines.append(f"  {col}=1 : {n_true:,}")
            lines.append(f"  {col}=0 : {n_false:,}")

        hbc_col = f"{mode}_has_bouquet_component"
        if hbc_col in result_df.columns:
            hbc = result_df[hbc_col].drop_nulls().cast(pl.Int64, strict=False)
            n_hbc = int((hbc == 1).sum())
            lines.append(f"  {hbc_col}=1 : {n_hbc:,}")

    return "\n".join(lines)



def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run WL-identifiability pipeline on a .smi file."
    )
    p.add_argument(
        "--smi",
        required=True,
        help="Path to the input .smi file (e.g. data/raw/ZINC/zinc250k.smi).",
    )
    p.add_argument(
        "--out",
        default=str(DEFAULT_OUT),
        help=f"Path to the shared SQLite database (default: {DEFAULT_OUT.relative_to(ROOT_DIR)}).",
    )
    p.add_argument(
        "--summary-csv",
        default=str(DEFAULT_SUMMARY_CSV),
        dest="summary_csv",
        help=f"Path to the dataset summary CSV (default: {DEFAULT_SUMMARY_CSV.relative_to(ROOT_DIR)}).",
    )
    p.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of parallel worker processes (default: 1). Use >1 for large datasets.",
    )
    p.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N molecules (default: all).",
    )
    p.add_argument(
        "--wl-steps",
        type=int,
        default=None,
        help="Fixed number of WL iterations. If omitted, estimated automatically from a sample.",
    )
    p.add_argument(
        "--sample-size",
        type=int,
        default=300,
        help="Number of molecules used for WL-step estimation (default: 300).",
    )
    p.add_argument(
        "--cap",
        type=int,
        default=50,
        help="Maximum WL iterations cap used during step estimation (default: 50).",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()

    smi_path = Path(args.smi)
    if not smi_path.exists():
        print(f"ERROR: file not found: {smi_path}")
        return 1

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    dataset_name = smi_path.name

    df = load_smi(smi_path, limit=args.limit)
    if df.is_empty():
        print("ERROR: no valid molecules loaded.")
        return 1

    k_max = 0
    k_p95 = 0

    if args.wl_steps is not None:
        fixed_wl_steps = args.wl_steps
        print(f"WL steps: {fixed_wl_steps} (user-specified)")
    else:
        print(f"Estimating WL steps from {min(args.sample_size, len(df)):,} molecules …")
        k_stats = estimate_fixed_wl_steps_from_dataframe(
            df, sample_size=args.sample_size, cap=args.cap
        )
        k_max = int(k_stats["K_max"])
        k_p95 = int(k_stats["K_p95"])
        fixed_wl_steps = k_p95
        print(f"WL steps: {fixed_wl_steps} (p95={k_p95}, max={k_max})")

    print(
        f"\nRunning pipeline on {len(df):,} molecules "
        f"(workers={args.workers}, wl_steps={fixed_wl_steps}) …\n"
    )
    t0 = time.perf_counter()
    result_df, n_failed = run_molecule_analysis_pipeline(
        df, fixed_wl_steps=fixed_wl_steps, n_workers=args.workers
    )
    elapsed = time.perf_counter() - t0

    conn = open_db(out_path)
    for row in result_df.to_dicts():
        try:
            insert_row(conn, row, dataset_name=dataset_name)
        except Exception as exc:
            print(f"[WARN] DB insert failed for {row.get('zinc_id')}: {exc}", file=sys.stderr)
    conn.close()

    summary_row = _build_summary_row(
        result_df=result_df,
        n_failed=n_failed,
        dataset_name=dataset_name,
        sample_size=args.sample_size,
        wl_cap=args.cap,
        k_max=k_max,
        k_p95=k_p95,
        wl_steps_used=fixed_wl_steps,
    )
    write_dataset_summary_row(args.summary_csv, summary_row)

    n_total = len(result_df)
    n_ok = n_total - n_failed
    summary_text = _build_summary_text(
        dataset_name=dataset_name,
        elapsed=elapsed,
        n_total=n_total,
        n_ok=n_ok,
        n_failed=n_failed,
        out_path=out_path,
        summary_csv=args.summary_csv,
        result_df=result_df,
    )
    print()
    print(summary_text)

    return 0 if n_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
