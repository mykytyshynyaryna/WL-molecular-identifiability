# Is 1-WL Expressivity Sufficient for Molecular Graphs?

A thesis project investigating whether the **1-Weisfeiler-Lehman (1-WL) graph isomorphism test** can uniquely identify every molecule in large chemical and biochemical graph datasets. The analysis is grounded in the theoretical framework of Kiefer (2020) — specifically the bouquet forest characterisation of non-identifiable graphs (Theorem 17).

---

## Research question

> Can the 1-WL algorithm distinguish every molecule in a standard benchmark dataset based solely on graph structure and atomic labels?

---

## Datasets

| Dataset | Domain | Graphs | Source |
|---------|--------|--------|--------|
| ZINC20 | Drug-like molecules (SMILES) | ~16 M | ZINC15 database |
| MUTAG | Mutagenic aromatic compounds | 188 | TU Dortmund |
| NCI1 | Anti-cancer activity screening | 4 110 | NCI DTP / TU Dortmund |
| NCI109 | Anti-cancer activity screening | 4 127 | NCI DTP / TU Dortmund |

---

## Key findings

Across **16,085,612 molecules** (MUTAG + NCI1 + NCI109 + ZINC20):

| Verdict | Count | Share |
|---------|-------|-------|
| Both WL modes accept (bouquet forest) | 15,990,031 | 99.41% |
| Atom-aware rejects, topological accepts | 62,988 | 0.39% |
| Both reject | 32,591 | 0.20% |
| Topological accepts, atom-aware rejects | 0 | 0.00% |

**Top rejection reasons:**

| Reason | Topological | Atom-aware |
|--------|-------------|------------|
| Component cycle ≠ C₅ | 58,242 | 28,997 |
| Petals not isomorphic | 34,089 | 2,416 |
| Multiple cycles in component | 3,248 | 1,178 |

See `notebooks/` for per-dataset breakdowns and visualisations of interesting cases.

---

## Repository structure

```
.
├── src/
│   └── wl_identifiability/        # Core algorithmic package
│       ├── wl.py                  # 1-WL algorithm
│       ├── graph_construction.py  # RDKit → NetworkX conversion
│       ├── flip_graph.py          # Flip graph construction
│       ├── skeleton.py            # Skeleton graph, Lemma 16 checks
│       ├── bouquet.py             # Bouquet detection (Theorem 17)
│       ├── experiments.py         # Pipeline orchestration
│       └── visualization.py       # WL coloring visualisation
├── examples/
│   └── run_pipeline.py            # Main CLI (ZINC-format .smi input)
├── scripts/
│   ├── database/db.py             # SQLite persistence layer
│   ├── pipeline/
│   │   ├── run_pipeline.py        # Richer CLI (--limit, --wl-steps, etc.)
│   │   └── compare_datasets.py    # Cross-dataset comparison
│   ├── download/
│   │   ├── download data/         # Per-dataset download scripts
│   │   └── parsing data/          # SMILES parsing scripts
│   ├── benchmarks/
│   │   └── benchmark_flip_graph.py
│   ├── profiling/
│   │   └── run all profiles/      # cProfile scripts + run_all_profiles.py
│   └── results/                   # Saved result summaries (markdown)
├── notebooks/
│   ├── interesting_cases_analysis_MUTAG.ipynb
│   ├── interesting_cases_analysis_NCI1.ipynb
│   ├── interesting_cases_analysis_NCI109.ipynb
│   └── interesting_cases_analysis_ZINC20.ipynb
├── tests/
│   ├── conftest.py
│   ├── test_wl.py
│   ├── test_flip_graph.py
│   ├── test_skeleton.py
│   ├── test_bouquet.py
│   ├── test_experiments.py
│   ├── test_visualization.py
│   └── test_graph_construction.py
├── data/                          # Dataset files (not tracked by git)
├── results/                       # Pipeline output (not tracked by git)
├── pyproject.toml
└── pixi.toml
```

---

## Setup

**Recommended — pixi (resolves rdkit via conda-forge automatically):**

```bash
pixi install
```

**Alternative — pip:**

```bash
pip install -e ".[dev]"
# If rdkit fails on pip:
conda install -c conda-forge rdkit
pip install -e ".[dev]" --no-deps
```

Python 3.10–3.12 required.

---

## Running the pipeline

### `examples/run_pipeline.py` — simple ZINC-format runner

```bash
python examples/run_pipeline.py --data data/AAAA.smi
python examples/run_pipeline.py --data data/AAAA.smi --workers 4 --sample 500
```

| Argument | Default | Description |
|---|---|---|
| `--data` | `data/AAAA.smi` | Path to `.smi` input file |
| `--sample` | `300` | Molecules sampled to estimate WL step count K |
| `--cap` | `50` | Max WL iterations during K estimation |
| `--db` | `results/all_results.db` | Shared SQLite DB (all datasets, one file) |
| `--summary-csv` | `results/summary.csv` | One summary row per dataset |
| `--workers` | `1` | Parallel worker processes |

### `scripts/pipeline/run_pipeline.py` — richer CLI

```bash
python scripts/pipeline/run_pipeline.py --smi data/ZINC20/zinc250k.smi --workers 4
python scripts/pipeline/run_pipeline.py --smi data/MUTAG.smi --limit 1000 --wl-steps 3
```

| Argument | Default | Description |
|---|---|---|
| `--smi` | required | Path to `.smi` file |
| `--out` | `results/all_results.db` | SQLite DB path |
| `--workers` | `1` | Parallel worker processes |
| `--limit` | all | Process only first N molecules |
| `--wl-steps` | auto | Fixed WL step count (skips estimation) |
| `--sample-size` | `300` | Sample size for WL step estimation |
| `--cap` | `50` | Max WL iterations cap |

**Run all datasets in sequence:**

```bash
for f in data/*.smi; do
    python examples/run_pipeline.py --data "$f"
done
```

---

## Input format

Space-separated `.smi` file with a header row:

```
smiles zinc_id
CO[C@H]1OC[C@@H](O)[C@H](O)[C@H]1O ZINC4371221
CC(=O)Nc1ccc(O)cc1 ZINC4509732
```

The `data/` directory is not tracked by git. Download scripts are in `scripts/download/`.

---

## Downloading datasets

```bash
# Download everything
python "scripts/download/download data/download_all.py"

# Individual datasets
python "scripts/download/download data/download_mutag.py"
python "scripts/download/download data/download_nci1.py"
python "scripts/download/download data/download_nci109.py"
python "scripts/download/download data/download_zinc.py"
```

See [scripts/download/README_downloads.md](scripts/download/README_downloads.md) for details.

---

## Tests

```bash
pixi run test
# or
pytest tests/ -v --tb=short
```

The suite covers WL coloring, graph construction, flip graph, skeleton, bouquet detection, experiments, visualisation (~100 tests across 8 files).

---

## Profiling

cProfile-based scripts in `scripts/profiling/run all profiles/`.

```bash
# Run all scenarios and print a comparison table
python "scripts/profiling/run all profiles/run_all_profiles.py"

# Individual scenarios
python "scripts/profiling/run all profiles/profile_bouquet_baseline.py"
python "scripts/profiling/run all profiles/profile_bouquet_optimized.py"
python "scripts/profiling/run all profiles/profile_flip_graph_baseline.py"
python "scripts/profiling/run all profiles/profile_flip_graph_optimized.py"
```

Common arguments: `--data data/AAAA.smi`, `--max-mols N`, `--reps N`, `--out profiles/`, `--no-save`.

See [scripts/profiling/README.md](scripts/profiling/README.md) for details.

---

## Benchmarks

```bash
python scripts/benchmarks/benchmark_flip_graph.py
python scripts/benchmarks/benchmark_flip_graph.py --data data/AAAC.smi --reps 200 --max-mols 500
```

---

## Notebooks

Jupyter notebooks in `notebooks/` contain per-dataset analysis of non-identifiable molecules (cases where the bouquet forest criterion fails or gives unexpected results).

```bash
jupyter notebook notebooks/
```

---

## Output

```
results/
├── all_results.db   # SQLite — one row per molecule, tagged by dataset_name
└── summary.csv      # One summary row per dataset (replaced on rerun)
```

Key columns in `all_results.db`:

| Column | Description |
|--------|-------------|
| `molecule_id` | Dataset identifier (e.g. ZINC ID) |
| `dataset_name` | Source file name |
| `smiles` | Input SMILES string |
| `n_nodes` / `n_edges` | Graph size |
| `n_colors_top` / `n_colors_atom` | WL color counts (topological / atom-aware) |
| `top_bouquet_forest_verdict` | `1` = bouquet forest (topological WL) |
| `atom_bouquet_forest_verdict` | `1` = bouquet forest (atom-aware WL) |

---

## Module overview

| Module | Responsibility |
|--------|---------------|
| `wl.py` | 1-WL coloring (fixed steps, converge, bounded modes) |
| `graph_construction.py` | RDKit `Mol` → NetworkX `Graph` with atom/bond attributes |
| `flip_graph.py` | Flip graph from WL partition via majority-edge rule |
| `skeleton.py` | Skeleton graph S_G; Notation 15 relation classifier; Lemma 16 validator |
| `bouquet.py` | C₅ detection, petal isomorphism, bouquet forest classification |
| `experiments.py` | `estimate_fixed_wl_steps_from_dataframe`, `analyze_single_molecule`, `run_molecule_analysis_pipeline` |
| `visualization.py` | WL coloring visualisation; `inspect_wl_behavior` |

---

## Notes

- **Package location:** `src/wl_identifiability/` (src layout — install before importing).
- **WL mode:** flip graph is built from atom-aware WL coloring, not topological.
- **Non-identifiability criterion:** at least two flip-graph components must be isomorphic bouquets (Definition 9 / Theorem 17, Kiefer 2020).
- **Multiprocessing:** `--workers > 1` uses `multiprocessing.Pool`; result order is non-deterministic.
- **Namespace conflict:** the root-level `wl_identifiability/` directory is a stale leftover — do not add `__init__.py` to it.

---

## Troubleshooting

**`ImportError: No module named 'wl_identifiability'`** — run `pixi install` or `pip install -e .`.

**`ImportError: No module named 'rdkit'`** — use `pixi install` or `conda install -c conda-forge rdkit`.

**`FileNotFoundError: data/AAAA.smi`** — place the file in `data/` or pass `--data`/`--smi` explicitly.
