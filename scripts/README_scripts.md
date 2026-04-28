# Scripts

All scripts run from the **project root** directory.

---

## Structure

```
scripts/
├── database/        SQLite persistence utilities
├── pipeline/        Full pipeline runners
├── download/        Dataset download and parsing scripts
├── benchmarks/      Timing and memory benchmarks
└── profiling/       cProfile-based profiling scripts
```

---

## database/

SQLite persistence used by the pipeline.

| File | Purpose |
|------|---------|
| `db.py` | Schema, `open_db`, `insert_row`, `write_dataset_summary_row` |

Not a standalone script — import from pipeline scripts:

```python
from scripts.database.db import open_db, insert_row
```

---

## pipeline/

Full pipeline runners.

| Script | Description |
|--------|-------------|
| `run_pipeline.py` | Richer CLI: `--limit`, `--wl-steps`, `--sample-size`, etc. |
| `compare_datasets.py` | Cross-dataset result comparison |

```bash
python scripts/pipeline/run_pipeline.py --smi data/AAAA.smi
python scripts/pipeline/run_pipeline.py --smi data/AAAA.smi --workers 4 --limit 10000
```

---

## download/

Download and parse benchmark datasets.

| Script | Dataset | Source |
|--------|---------|--------|
| `download data/download_all.py` | All datasets | runs all scripts in sequence |
| `download data/download_mutag.py` | MUTAG | TU Dortmund |
| `download data/download_nci1.py` | NCI1 | TU Dortmund |
| `download data/download_nci109.py` | NCI109 | TU Dortmund |
| `download data/download_zinc.py` | ZINC 250k | GitHub |
| `download data/download_zinc_from_uri.py` | ZINC (custom tranche) | ZINC15 URI list |
| `parsing data/parse_mutag_to_smiles.py` | MUTAG → SMILES | — |
| `parsing data/parse_nci1_to_smiles.py` | NCI1 → SMILES | — |
| `parsing data/parse_nci109_to_smiles.py` | NCI109 → SMILES | — |

```bash
python "scripts/download/download data/download_all.py"
python "scripts/download/download data/download_mutag.py"
python "scripts/download/download data/download_zinc_from_uri.py" --uri path/to/file.uri
```

All files are saved to `data/raw/<DATASET_NAME>/`. See [download/README_downloads.md](download/README_downloads.md) for details.

---

## benchmarks/

Head-to-head timing and memory benchmarks.

| Script | What it benchmarks |
|--------|-------------------|
| `benchmark_flip_graph.py` | Baseline (NumPy matrix) vs optimised (edge-list) flip graph construction |

```bash
python scripts/benchmarks/benchmark_flip_graph.py
python scripts/benchmarks/benchmark_flip_graph.py --data data/AAAC.smi --reps 200 --max-mols 500
python scripts/benchmarks/benchmark_flip_graph.py --no-save
```

---

## profiling/

cProfile-based profiling of algorithm implementations.

| Script | What it profiles |
|--------|-----------------|
| `profile_bouquet_baseline.py` | Bouquet detection — baseline (`nx.is_isomorphic`) |
| `profile_bouquet_optimized.py` | Bouquet detection — optimised (AHU signatures) |
| `profile_flip_graph_baseline.py` | Flip graph — dense NumPy approach |
| `profile_flip_graph_optimized.py` | Flip graph — edge-list approach |
| `run_all_profiles.py` | Runs all four and prints a comparison table |

```bash
python "scripts/profiling/run all profiles/run_all_profiles.py"
python "scripts/profiling/run all profiles/profile_bouquet_optimized.py" --data data/AAAC.smi --reps 10
```

See [profiling/README.md](profiling/README.md) for full details.
