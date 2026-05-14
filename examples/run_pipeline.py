from __future__ import annotations

import argparse
import sqlite3
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import polars as pl

from scripts.database.db import (
    insert_row as db_insert_row,
)
from scripts.database.db import (
    open_db,
    write_dataset_summary_row,
)
from wl_identifiability.experiments import (
    estimate_fixed_wl_steps_from_dataframe,
    run_molecule_analysis_pipeline,
)

DEFAULT_RESULTS_DIR = Path("results")
DEFAULT_GLOBAL_DB = DEFAULT_RESULTS_DIR / "all_results.db"
DEFAULT_SUMMARY_CSV = DEFAULT_RESULTS_DIR / "summary.csv"


def load_dataset(path: str) -> pl.DataFrame:
    df = pl.read_csv(path, separator=" ", has_header=True)
    if "zinc_id" in df.columns:
        df = df.select(["smiles", "zinc_id"]).rename({"zinc_id": "molecule_id"})
    elif "id" in df.columns:
        df = df.select(["smiles", "id"]).rename({"id": "molecule_id"})
    else:
        raise ValueError(f"No 'id' or 'zinc_id' column found in {path}")
    df = df.with_columns(pl.col("molecule_id").cast(pl.Utf8))
    print(f"  Loaded {len(df):,} rows from {path}")
    return df


def check_imports() -> bool:
    modules = [
        "wl_identifiability",
        "wl_identifiability.wl",
        "wl_identifiability.graph_construction",
        "wl_identifiability.flip_graph",
        "wl_identifiability.bouquet",
        "wl_identifiability.visualization",
        "wl_identifiability.experiments",
    ]
    ok = True
    for name in modules:
        try:
            __import__(name)
            print(f"  [OK]   {name}")
        except Exception:
            print(f"  [FAIL] {name}")
            traceback.print_exc()
            ok = False
    return ok


def _safe_int(value: object) -> int:
    return 0 if value is None else int(value)


def _safe_float(value: object) -> float:
    return 0.0 if value is None else float(value)


def build_summary_text(
    res: pl.DataFrame,
    bad: int,
    data_path: Path,
    sample_size: int,
    cap: int,
    k_max: int,
    k_p95: int,
) -> str:
    total = len(res)
    ok_count = max(total - bad, 0)
    ok_pct = (100.0 * ok_count / total) if total > 0 else 0.0

    def _col_sum(col: str) -> int:
        if col not in res.columns:
            return 0
        return _safe_int(res.select(pl.col(col).cast(pl.Int64, strict=False).sum()).item())

    def _col_mean(col: str) -> float:
        if col not in res.columns:
            return 0.0
        return _safe_float(res.select(pl.col(col).mean()).item())

    def _col_max(col: str) -> int:
        if col not in res.columns:
            return 0
        v = res.select(pl.col(col).max()).item()
        return _safe_int(v)

    top_bouquet_forest = _col_sum("top_bouquet_forest_verdict")
    atom_bouquet_forest = _col_sum("atom_bouquet_forest_verdict")
    avg_colors_top = _col_mean("n_colors_top")
    avg_colors_atom = _col_mean("n_colors_atom")
    avg_color_ratio = _col_mean("color_ratio_atom_to_top")
    avg_flipped_edges = _col_mean("n_flipped_edges")
    max_flipped_edges = _col_max("n_flipped_edges")

    lines = [
        "=== Pipeline summary ===",
        f"Dataset           : {data_path}",
        f"Sample size       : {sample_size}",
        f"WL cap            : {cap}",
        f"K_max             : {k_max}",
        f"K_p95             : {k_p95}",
        "",
        f"Total molecules   : {total:>8,}",
        f"Parsed OK         : {ok_count:>8,}  ({ok_pct:.1f}%)",
        f"Parse/WL errors   : {bad:>8,}",
        f"BF verdict (top)  : {top_bouquet_forest:>8,}",
        f"BF verdict (atom) : {atom_bouquet_forest:>8,}",
        f"Avg colors (top)  : {avg_colors_top:>11.3f}",
        f"Avg colors (atom) : {avg_colors_atom:>11.3f}",
        f"Avg color ratio   : {avg_color_ratio:>11.3f}",
        f"Avg flipped edges : {avg_flipped_edges:>11.3f}",
        f"Max flipped edges : {max_flipped_edges:>8,}",
    ]

    return "\n".join(lines)


def build_summary_row(
    res: pl.DataFrame,
    bad: int,
    dataset_name: str,
    sample_size: int,
    cap: int,
    k_max: int,
    k_p95: int,
) -> dict:
    """Build one structured summary row for dataset_summaries.csv.

    Field names and order must stay in sync with SUMMARY_FIELDS in
    scripts/database/db.py so the CSV header is stable across datasets.
    """
    total = len(res)
    ok_count = max(total - bad, 0)
    ok_pct = round(100.0 * ok_count / total, 1) if total > 0 else 0.0

    def _col_sum(col: str) -> int:
        if col not in res.columns:
            return 0
        return _safe_int(res.select(pl.col(col).cast(pl.Int64, strict=False).sum()).item())

    def _col_mean(col: str) -> float:
        if col not in res.columns:
            return 0.0
        return round(_safe_float(res.select(pl.col(col).mean()).item()), 3)

    def _col_max(col: str) -> int:
        if col not in res.columns:
            return 0
        return _safe_int(res.select(pl.col(col).max()).item())

    return {
        "dataset_name": dataset_name,
        "sample_size": sample_size,
        "wl_cap": cap,
        "k_max": k_max,
        "k_p95": k_p95,
        "wl_steps_used": k_max,
        "total_molecules": total,
        "parsed_ok": ok_count,
        "parsed_ok_pct": ok_pct,
        "parse_wl_errors": bad,
        "bf_verdict_top": _col_sum("top_bouquet_forest_verdict"),
        "bf_verdict_atom": _col_sum("atom_bouquet_forest_verdict"),
        "avg_colors_top": _col_mean("n_colors_top"),
        "avg_colors_atom": _col_mean("n_colors_atom"),
        "avg_color_ratio": _col_mean("color_ratio_atom_to_top"),
        "avg_flipped_edges": _col_mean("n_flipped_edges"),
        "max_flipped_edges": _col_max("n_flipped_edges"),
    }


def print_summary(
    res: pl.DataFrame,
    bad: int,
    data_path: Path,
    sample_size: int,
    cap: int,
    k_max: int,
    k_p95: int,
) -> str:
    summary_text = build_summary_text(
        res=res,
        bad=bad,
        data_path=data_path,
        sample_size=sample_size,
        cap=cap,
        k_max=k_max,
        k_p95=k_p95,
    )
    print("\n" + summary_text)
    return summary_text


def finalize_sqlite_db(db_path: Path) -> None:
    """
    Force SQLite to leave WAL mode so that .db-wal and .db-shm disappear.
    Safe to call after the pipeline has finished writing.
    """
    if not db_path.exists():
        return

    with sqlite3.connect(str(db_path)) as conn:
        conn.execute("PRAGMA wal_checkpoint(FULL);")
        conn.execute("PRAGMA journal_mode=DELETE;")
        conn.commit()

    wal_path = db_path.with_name(db_path.name + "-wal")
    shm_path = db_path.with_name(db_path.name + "-shm")

    for extra_file in (wal_path, shm_path):
        try:
            if extra_file.exists():
                extra_file.unlink()
        except PermissionError:
            print(f"  [WARN] Could not remove {extra_file.name} — file is still open in another program.")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run the full molecule identifiability pipeline.")
    p.add_argument("--data", default="data/AAAA.smi", help="Path to .smi dataset file")
    p.add_argument("--sample", type=int, default=300, help="Sample size for WL step estimation")
    p.add_argument("--cap", type=int, default=50, help="Max WL iterations cap")
    p.add_argument(
        "--db",
        default=None,
        help=(
            "SQLite DB path for per-molecule results. "
            "Default: results/<dataset_stem>.db (derived from --data filename)."
        ),
    )
    p.add_argument(
        "--summary-csv",
        dest="summary_csv",
        default=str(DEFAULT_SUMMARY_CSV),
        help=(
            "Shared CSV with one summary row per dataset. "
            f"Default: {DEFAULT_SUMMARY_CSV}. Rows for a re-run dataset are replaced."
        ),
    )
    p.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of worker processes (default: 1 = single-threaded)",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()

    print("=== 1. Checking module imports ===")
    if not check_imports():
        print("\nAborting: one or more modules failed to import.")
        return 1

    print("\n=== 2. Loading dataset ===")
    data_path = Path(args.data)
    if not data_path.exists():
        print(f"  ERROR: file not found: {data_path}")
        return 1

    df = load_dataset(str(data_path))
    if len(df) == 0:
        print("  ERROR: dataset is empty after parsing.")
        return 1

    dataset_name = data_path.name

    db_path = Path(args.db) if args.db else DEFAULT_RESULTS_DIR / f"{data_path.stem}.db"
    db_path.parent.mkdir(parents=True, exist_ok=True)

    summary_csv_path = Path(args.summary_csv)
    summary_csv_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"  SQLite DB        : {db_path}")
    print(f"  Summary CSV      : {summary_csv_path}  (one row per dataset)")

    print(f"\n=== 3. Estimating WL steps (sample={args.sample}, cap={args.cap}) ===")
    k_stats = estimate_fixed_wl_steps_from_dataframe(
        df,
        sample_size=args.sample,
        cap=args.cap,
    )
    K = int(k_stats["K_max"])
    K_p95 = int(k_stats["K_p95"])

    print(f"  K_max = {K}")
    print(f"  K_p95 = {K_p95}")

    print(f"\n=== 4. Running pipeline on {len(df):,} molecules (K={K}, workers={args.workers}) ===")
    res, bad = run_molecule_analysis_pipeline(
        df,
        fixed_wl_steps=K,
        n_workers=args.workers,
    )

    print(f"\n=== 5. Writing results to SQLite ({db_path}) ===")
    conn = open_db(db_path)
    for row in res.iter_rows(named=True):
        db_insert_row(conn, row, dataset_name=dataset_name)
    conn.close()

    finalize_sqlite_db(db_path)

    print_summary(
        res=res,
        bad=bad,
        data_path=data_path,
        sample_size=args.sample,
        cap=args.cap,
        k_max=K,
        k_p95=K_p95,
    )

    summary_row = build_summary_row(
        res=res,
        bad=bad,
        dataset_name=dataset_name,
        sample_size=args.sample,
        cap=args.cap,
        k_max=K,
        k_p95=K_p95,
    )
    write_dataset_summary_row(summary_csv_path, summary_row)

    print(f"\n  Summary CSV row  : {summary_csv_path}  (dataset_name={dataset_name})")
    print(f"  SQLite DB saved  : {db_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
