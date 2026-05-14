"""Minimal SQLite persistence layer.

Only the master process calls open_db / insert_row.
Workers return plain dicts; no DB imports needed there.

Schema migration strategy
--------------------------
SQLite's ``CREATE TABLE IF NOT EXISTS`` is a no-op when the table already
exists, so adding new columns to SCHEMA does *not* update old .db files.
``_migrate_results_schema`` closes this gap: it reads the live column list
via PRAGMA and issues ``ALTER TABLE … ADD COLUMN`` for any column that is
present in EXPECTED_RESULTS_COLUMNS but absent from the table.  Adding a
nullable column is always safe and non-destructive — old rows get NULL for
the new fields, which is the correct sentinel for "not yet computed".

Old columns that may exist in legacy files (bouquet_forest_verdict, reason,
mode, top_non_identifiable, atom_non_identifiable) are left untouched —
SQLite cannot drop columns portably, and keeping them is harmless.
"""

from __future__ import annotations

import csv
import json
import sqlite3
from pathlib import Path

SCHEMA = """
CREATE TABLE IF NOT EXISTS results (
    molecule_id                  TEXT    NOT NULL PRIMARY KEY,
    dataset_name                 TEXT,   -- source dataset file name (e.g. AAAC.smi)
    smiles                       TEXT    NOT NULL,
    n_nodes                      INTEGER,
    n_edges                      INTEGER,
    wl_stats                     TEXT,   -- JSON: wl_converged/iters/n_colors for top and atom
    flip_graph_stats             TEXT,   -- JSON: within/between copy/flip counts (atom flip graph)
    top_bouquet_forest_verdict   INTEGER,-- 0 / 1 / NULL  (topology-only WL flip graph)
    atom_bouquet_forest_verdict  INTEGER,-- 0 / 1 / NULL  (atom-aware WL flip graph)
    top_bouquet_forest_reason    TEXT,   -- short reason code for topology verdict
    atom_bouquet_forest_reason   TEXT,   -- short reason code for atom verdict
    top_has_bouquet_component    INTEGER,-- 0 / 1 / NULL  1 if ≥1 component of topology flip graph is a bouquet
    atom_has_bouquet_component   INTEGER -- 0 / 1 / NULL  1 if ≥1 component of atom-aware flip graph is a bouquet
);
"""

EXPECTED_RESULTS_COLUMNS: list[tuple[str, str]] = [
    ("molecule_id", "TEXT"),
    ("dataset_name", "TEXT"),
    ("smiles", "TEXT"),
    ("n_nodes", "INTEGER"),
    ("n_edges", "INTEGER"),
    ("wl_stats", "TEXT"),
    ("flip_graph_stats", "TEXT"),
    ("top_bouquet_forest_verdict", "INTEGER"),
    ("atom_bouquet_forest_verdict", "INTEGER"),
    ("top_bouquet_forest_reason", "TEXT"),
    ("atom_bouquet_forest_reason", "TEXT"),
    ("top_has_bouquet_component", "INTEGER"),
    ("atom_has_bouquet_component", "INTEGER"),
]


SUMMARY_FIELDS: list[str] = [
    "dataset_name",
    "sample_size",
    "wl_cap",
    "k_max",
    "k_p95",
    "wl_steps_used",
    "total_molecules",
    "parsed_ok",
    "parsed_ok_pct",
    "parse_wl_errors",
    "bf_verdict_top",
    "bf_verdict_atom",
    "avg_colors_top",
    "avg_colors_atom",
    "avg_color_ratio",
    "avg_flipped_edges",
    "max_flipped_edges",
]

SCHEMA_GRAPHS = """
CREATE TABLE IF NOT EXISTS graph_results (
    graph_id                INTEGER NOT NULL PRIMARY KEY,
    graph_label             INTEGER,
    dataset                 TEXT,
    n_nodes                 INTEGER,
    n_edges                 INTEGER,
    wl_iterations           INTEGER,
    wl_converged            INTEGER,-- 0 / 1 / NULL
    n_colors_top            INTEGER,
    flip_nodes              INTEGER,
    flip_edges              INTEGER,
    n_flipped_edges         INTEGER,
    is_bouquet_forest       INTEGER,-- 0 / 1 / NULL
    non_identifiable        INTEGER,-- 0 / 1 / NULL
    num_bouquets            INTEGER,
    reason                  TEXT,
    error                   TEXT
);
"""


def _migrate_results_schema(conn: sqlite3.Connection) -> None:
    """
    Ensure the `results` table has every column in EXPECTED_RESULTS_COLUMNS.

    For each missing column, issues ``ALTER TABLE results ADD COLUMN …``.
    Adding a nullable column is always safe: existing rows receive NULL,
    which is the correct default for "not yet computed".

    Old columns that are no longer in the expected list (e.g. the former
    ``bouquet_forest_verdict``, ``reason``, ``mode``) are left in place —
    SQLite cannot drop columns portably, and keeping them is harmless.

    Raises RuntimeError with a clear message if the table exists but is
    missing a NOT NULL column without a default (i.e. PRIMARY KEY columns),
    because those cannot be added via ALTER TABLE.
    """
    rows = conn.execute("PRAGMA table_info(results);").fetchall()
    existing = {row[1] for row in rows}

    non_nullable = {"molecule_id", "smiles"}

    for col_name, col_type in EXPECTED_RESULTS_COLUMNS:
        if col_name in existing:
            continue
        if col_name in non_nullable:
            raise RuntimeError(
                f"Database schema is incompatible: column '{col_name}' is missing "
                f"from the 'results' table and cannot be added automatically because "
                f"it is NOT NULL.  Please delete or rename the existing .db file so "
                f"a fresh database is created with the correct schema."
            )
        conn.execute(f"ALTER TABLE results ADD COLUMN {col_name} {col_type};")
        print(f"  [db] migrated: added column '{col_name}' ({col_type}) to results table.")


def open_db(path: str | Path) -> sqlite3.Connection:
    """
    Open (or create) the SQLite DB and ensure the molecule results schema is current.

    On a brand-new file the full schema is created from SCHEMA.
    On an existing file ``_migrate_results_schema`` adds any columns that were
    introduced since the file was first created (safe ALTER TABLE ADD COLUMN).
    """
    conn = sqlite3.connect(str(path), isolation_level=None)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute(SCHEMA)
    _migrate_results_schema(conn)
    return conn


def open_graph_db(path: str | Path) -> sqlite3.Connection:
    """Open (or create) the SQLite DB and ensure the graph-dataset schema exists."""
    conn = sqlite3.connect(str(path), isolation_level=None)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute(SCHEMA_GRAPHS)
    return conn


def insert_graph_row(conn: sqlite3.Connection, row: dict, dataset: str = "") -> None:
    """Insert one result dict produced by analyze_one_graph in run_reddit_pipeline."""
    verdict = row.get("is_bouquet_forest")
    non_id = row.get("non_identifiable")
    wl_conv = row.get("wl_converged")

    conn.execute(
        """
        INSERT OR REPLACE INTO graph_results
            (graph_id, graph_label, dataset,
             n_nodes, n_edges,
             wl_iterations, wl_converged, n_colors_top,
             flip_nodes, flip_edges, n_flipped_edges,
             is_bouquet_forest, non_identifiable, num_bouquets,
             reason, error)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            row.get("graph_id"),
            row.get("graph_label"),
            dataset,
            row.get("n_nodes"),
            row.get("n_edges"),
            row.get("wl_iterations"),
            int(wl_conv) if wl_conv is not None else None,
            row.get("n_colors_top"),
            row.get("flip_nodes"),
            row.get("flip_edges"),
            row.get("n_flipped_edges"),
            int(verdict) if verdict is not None else None,
            int(non_id) if non_id is not None else None,
            row.get("num_bouquets"),
            row.get("reason"),
            row.get("error"),
        ),
    )


def insert_row(conn: sqlite3.Connection, row: dict, dataset_name: str = "") -> None:
    """Insert one result dict produced by analyze_single_molecule.

    ``dataset_name`` should be the source file name (e.g. ``AAAC.smi``).
    It defaults to an empty string for backward compatibility with callers
    that do not pass this argument.
    """
    wl_keys = {
        "wl_converged_top",
        "wl_iters_top",
        "wl_budget_top",
        "n_colors_top",
        "wl_converged_atom",
        "wl_iters_atom",
        "wl_budget_atom",
        "n_colors_atom",
        "color_ratio_atom_to_top",
    }
    flip_keys = {
        "within_copy",
        "within_flip",
        "between_copy",
        "between_flip",
        "n_flipped_edges",
    }

    wl_stats = {k: row[k] for k in wl_keys if k in row}
    flip_stats = {k: row[k] for k in flip_keys if k in row}

    def _int_or_none(v):
        return int(v) if v is not None else None

    conn.execute(
        """
        INSERT OR REPLACE INTO results
            (molecule_id, dataset_name, smiles, n_nodes, n_edges,
             wl_stats, flip_graph_stats,
             top_bouquet_forest_verdict,  atom_bouquet_forest_verdict,
             top_bouquet_forest_reason,   atom_bouquet_forest_reason,
             top_has_bouquet_component,   atom_has_bouquet_component)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            str(row.get("molecule_id", "")),
            dataset_name,
            str(row.get("smiles", "")),
            row.get("n_nodes"),
            row.get("n_edges"),
            json.dumps(wl_stats) if wl_stats else None,
            json.dumps(flip_stats) if flip_stats else None,
            _int_or_none(row.get("top_bouquet_forest_verdict")),
            _int_or_none(row.get("atom_bouquet_forest_verdict")),
            row.get("top_bouquet_forest_reason"),
            row.get("atom_bouquet_forest_reason"),
            _int_or_none(row.get("top_has_bouquet_component")),
            _int_or_none(row.get("atom_has_bouquet_component")),
        ),
    )


def write_dataset_summary_row(csv_path: str | Path, row: dict) -> None:
    """Append (or replace) a dataset summary row in the CSV summary file.

    If the file does not exist it is created with a header.
    If a row for the same ``dataset_name`` already exists it is replaced so
    that re-running a dataset does not produce duplicate entries.

    ``row`` must contain at least the keys defined in ``SUMMARY_FIELDS``.
    Extra keys are silently ignored; missing keys are written as empty strings.
    """
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    dataset_name = str(row.get("dataset_name", ""))

    existing_rows: list[dict] = []
    if csv_path.exists():
        with csv_path.open(newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            for existing_row in reader:
                if existing_row.get("dataset_name") != dataset_name:
                    existing_rows.append(existing_row)

    new_row = {field: row.get(field, "") for field in SUMMARY_FIELDS}
    existing_rows.append(new_row)

    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(existing_rows)
