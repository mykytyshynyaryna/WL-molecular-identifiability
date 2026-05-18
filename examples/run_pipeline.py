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
    # Force SQLite out of WAL mode so .db-wal/.db-shm files are cleaned up.
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


def run_on_file(
    data_path: Path,
    db_path: Path | None,
    summary_csv_path: Path,
    sample: int,
    cap: int,
    workers: int,
    indices: list[int] | None,
) -> int:
    if not data_path.exists():
        print(f"  ERROR: file not found: {data_path}")
        return 1

    df = load_dataset(str(data_path))
    if len(df) == 0:
        print("  ERROR: dataset is empty after parsing.")
        return 1

    if indices is not None:
        valid = [i for i in indices if 0 <= i < len(df)]
        skipped = len(indices) - len(valid)
        if skipped:
            print(f"  [WARN] {skipped} index/indices out of range (dataset has {len(df)} rows), skipped")
        df = df[valid]
        if len(df) == 0:
            print("  ERROR: no rows remain after applying --indices filter.")
            return 1
        print(f"  Filtered to {len(df):,} rows via --indices")

    dataset_name = data_path.name
    resolved_db = db_path if db_path else DEFAULT_RESULTS_DIR / f"{data_path.stem}.db"
    resolved_db.parent.mkdir(parents=True, exist_ok=True)
    summary_csv_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"  SQLite DB        : {resolved_db}")
    print(f"  Summary CSV      : {summary_csv_path}  (one row per dataset)")

    print(f"\n  Estimating WL steps (sample={sample}, cap={cap}) ...")
    k_stats = estimate_fixed_wl_steps_from_dataframe(df, sample_size=sample, cap=cap)
    K = int(k_stats["K_max"])
    K_p95 = int(k_stats["K_p95"])
    print(f"  K_max = {K}")
    print(f"  K_p95 = {K_p95}")

    print(f"\n  Running pipeline on {len(df):,} molecules (K={K}, workers={workers}) ...")
    res, bad = run_molecule_analysis_pipeline(df, fixed_wl_steps=K, n_workers=workers)

    print(f"\n  Writing results to SQLite ({resolved_db}) ...")
    conn = open_db(resolved_db)
    for row in res.iter_rows(named=True):
        db_insert_row(conn, row, dataset_name=dataset_name)
    conn.close()

    finalize_sqlite_db(resolved_db)

    print_summary(
        res=res,
        bad=bad,
        data_path=data_path,
        sample_size=sample,
        cap=cap,
        k_max=K,
        k_p95=K_p95,
    )

    summary_row = build_summary_row(
        res=res,
        bad=bad,
        dataset_name=dataset_name,
        sample_size=sample,
        cap=cap,
        k_max=K,
        k_p95=K_p95,
    )
    write_dataset_summary_row(summary_csv_path, summary_row)

    print(f"\n  Summary CSV row  : {summary_csv_path}  (dataset_name={dataset_name})")
    print(f"  SQLite DB saved  : {resolved_db}")
    return 0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run the full molecule identifiability pipeline.")
    input_group = p.add_mutually_exclusive_group()
    input_group.add_argument("--data", default=None, help="Path to a single .smi dataset file")
    input_group.add_argument(
        "--data-dir",
        dest="data_dir",
        default=None,
        help="Directory containing .smi files; all are processed in sorted order",
    )
    p.add_argument(
        "--indices",
        default=None,
        help="Comma-separated 0-based row indices to process (e.g. --indices 0,5,10); only valid with --data",
    )
    p.add_argument("--sample", type=int, default=300, help="Sample size for WL step estimation")
    p.add_argument("--cap", type=int, default=50, help="Max WL iterations cap")
    p.add_argument(
        "--db",
        default=None,
        help=(
            "SQLite DB path for per-molecule results. "
            "Default: results/<dataset_stem>.db (derived from --data filename). "
            "Not valid with --data-dir."
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

    if args.data is None and args.data_dir is None:
        print("ERROR: provide --data <file.smi> or --data-dir <directory>", file=sys.stderr)
        return 1

    if args.data_dir and args.db:
        print("ERROR: --db cannot be used with --data-dir (each file gets its own database)", file=sys.stderr)
        return 1

    if args.data_dir and args.indices:
        print("ERROR: --indices cannot be used with --data-dir", file=sys.stderr)
        return 1

    print("=== 1. Checking module imports ===")
    if not check_imports():
        print("\nAborting: one or more modules failed to import.")
        return 1

    indices: list[int] | None = None
    if args.indices:
        try:
            indices = [int(x.strip()) for x in args.indices.split(",")]
        except ValueError:
            print("ERROR: --indices must be comma-separated integers, e.g. --indices 0,5,10", file=sys.stderr)
            return 1

    summary_csv_path = Path(args.summary_csv)

    if args.data_dir:
        data_dir = Path(args.data_dir)
        if not data_dir.is_dir():
            print(f"ERROR: --data-dir path is not a directory: {data_dir}", file=sys.stderr)
            return 1
        smi_files = sorted(data_dir.glob("*.smi"))
        if not smi_files:
            print(f"ERROR: no .smi files found in {data_dir}", file=sys.stderr)
            return 1
        print(f"\nFound {len(smi_files)} .smi file(s) in {data_dir}")
        failed: list[Path] = []
        for i, smi_path in enumerate(smi_files, 1):
            print(f"\n[{i}/{len(smi_files)}] Processing: {smi_path.name}")
            rc = run_on_file(
                data_path=smi_path,
                db_path=None,
                summary_csv_path=summary_csv_path,
                sample=args.sample,
                cap=args.cap,
                workers=args.workers,
                indices=None,
            )
            if rc != 0:
                failed.append(smi_path)
        if failed:
            print(f"\n{len(failed)} file(s) failed:")
            for f in failed:
                print(f"  {f.name}")
            return 1
        print(f"\nAll {len(smi_files)} file(s) processed successfully.")
        return 0

    return run_on_file(
        data_path=Path(args.data),
        db_path=Path(args.db) if args.db else None,
        summary_csv_path=summary_csv_path,
        sample=args.sample,
        cap=args.cap,
        workers=args.workers,
        indices=indices,
    )


if __name__ == "__main__":
    sys.exit(main())
