# Is 1-WL Expressivity Sufficient for Molecular Graphs?

> **Note:**
> - branch `main` — clean implementation branch with source code, examples, tests, and minimal reproducibility instructions.
> - branch `thesis` — extended research branch with profiling, benchmarks, notebooks, interesting cases, figures, and thesis-related experimental materials.

This project investigates whether the **1-Weisfeiler-Lehman (1-WL) graph isomorphism test** can uniquely identify every molecule in standard benchmark datasets. The pipeline runs WL coloring on molecular graphs, builds a flip graph from the stable color partition, and checks each connected component against the bouquet forest criterion to detect non-identifiable molecules.

---

## Background

**1-Weisfeiler-Lehman (1-WL) test** is an iterative colour-refinement algorithm on graphs. Each node starts with an initial colour (label); at every step every node aggregates the multiset of its neighbours' colours and hashes it into a new colour. The process stops when no colour changes. Two graphs are *non-isomorphic* by 1-WL if they produce different final colour histograms, and *potentially isomorphic* otherwise.

**Two WL modes used in this project:**
- *Topological* — initial node colour is uniform (structure only, no atom labels).
- *Atom-aware* — initial node colour encodes the atomic symbol, so chemical identity is used from the first iteration.

**Bouquet forest** (Kiefer 2020, Definition 9 / Theorem 17): a graph is a *bouquet* if it consists of a single C₅ cycle with exactly five rooted trees ("petals") — one attached at each cycle vertex — where all five petals are mutually isomorphic. A *bouquet forest* is a disjoint union of isomorphic bouquets. Theorem 17 states that a graph is *not* uniquely identified by 1-WL if and only if its flip graph decomposes into a bouquet forest with at least two components.

**Reference:** Kiefer, S. (2020). *The Weisfeiler-Leman Algorithm: Its Power and Limitations.* Habilitation thesis, RWTH Aachen University.

---

## Datasets

| Dataset | Domain | Graphs | Source |
|---------|--------|--------|--------|
| MUTAG | Mutagenic aromatic compounds | 188 | TU Dortmund |
| NCI1 | Anti-cancer activity screening | 4 110 | NCI DTP / TU Dortmund |
| NCI109 | Anti-cancer activity screening | 4 127 | NCI DTP / TU Dortmund |
| ZINC20 | Drug-like molecules (SMILES) | ~16 M | ZINC15 database |

---

## Key findings

Across **16,085,610 molecules** (MUTAG + NCI1 + NCI109 + ZINC20):

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

**Interpreting verdicts:** A molecule is called *identifiable* (verdict = 1) when its flip graph is a bouquet forest — meaning 1-WL is sufficient to distinguish it from all structurally similar molecules. "Both WL modes accept" means identifiable under both topological and atom-aware coloring. "Both reject" means 1-WL cannot distinguish the molecule regardless of whether atom labels are used.

---

## Prerequisites

- Python 3.12
- [pixi](https://pixi.sh) (recommended — handles rdkit automatically via conda-forge)

---

## 1. Install

**Recommended:**

```console
pixi install
```

This creates an isolated environment with all dependencies (networkx, numpy, polars, rdkit, matplotlib, pytest).

**Key dependency versions (resolved by pixi):**

| Package | Version |
|---------|---------|
| Python | 3.12 |
| RDKit | ≥ 2023.09 |
| NetworkX | ≥ 3.1 |
| polars | ≥ 1.0 |
| pytest | ≥ 7.0 |

**Alternative (pip):**

```console
pip install -e ".[dev]"
```

> If rdkit fails on pip, install it via conda first:
> ```console
> conda install -c conda-forge rdkit
> pip install -e ".[dev]" --no-deps
> ```

**Recommended first steps after install:**

```console
# 1. Verify installation
pixi run test

# 2. Quick smoke test
python examples/run_pipeline.py --data data/processed/MUTAG/mutag_smiles.smi --sample 50
```

---

## 2. Quick checks

### From the terminal (no Python script needed)

After installation the `wl-check` command is available inside the pixi environment:

```console
pixi run wl-check "CCCN"
pixi run wl-check "CCCN" "c1ccccc1" "CC(=O)O"
```

Output:

```
CCCN: identifiable
c1ccccc1: NOT identifiable
CC(=O)O: identifiable
```

Pass any number of SMILES strings as positional arguments. Exit code is `0` if all succeeded, `1` if any SMILES could not be parsed.

---

### As a Python library

After installation, the package exposes three simple functions — no files, no CLI needed:

```python
from wl_identifiability import is_smi_identifiable, is_mol_identifiable, is_graph_identifiable
```

### Check a single SMILES string

```python
from wl_identifiability import is_smi_identifiable

print(is_smi_identifiable("CCCN"))          # True — identifiable
print(is_smi_identifiable("C1=CC=CC=C1"))   # True — benzene, identifiable
```

Returns `True` if 1-WL (atom-aware) is sufficient to uniquely identify the molecule, `False` otherwise. Raises `ValueError` if RDKit cannot parse the SMILES.

### Check an RDKit molecule object

```python
from rdkit import Chem
from wl_identifiability import is_mol_identifiable

mol = Chem.MolFromSmiles("CCCN")
print(is_mol_identifiable(mol))   # True
```

### Check a NetworkX graph (topological mode)

```python
import networkx as nx
from wl_identifiability import is_graph_identifiable

G = nx.cycle_graph(6)
print(is_graph_identifiable(G))   # True
```

Uses topological 1-WL (no atom labels) — suitable for any NetworkX graph.

### Check a list of molecules

```python
from wl_identifiability import is_smi_identifiable

smiles_list = ["CCCN", "c1ccccc1", "CC(=O)O"]
results = {smi: is_smi_identifiable(smi) for smi in smiles_list}
print(results)
# {'CCCN': True, 'c1ccccc1': True, 'CC(=O)O': True}
```

---

## 3. Download datasets

> **Ready to test right away:** `data/processed/MUTAG/mutag_smiles.smi` (188 molecules, 8 KB) is already included in the repository — no download or parsing needed. Jump straight to [step 5](#5-run-the-pipeline):
> ```console
> python examples/run_pipeline.py --data data/processed/MUTAG/mutag_smiles.smi
> ```

All full datasets go into `data/raw/`. Run from the project root.

**Download everything at once:**

```console
python scripts/download/download_data/download_all.py
```

**Or download individually:**

```console
python scripts/download/download_data/download_mutag.py
python scripts/download/download_data/download_nci1.py
python scripts/download/download_data/download_nci109.py
python scripts/download/download_data/download_zinc.py
```

After downloading, `data/raw/` will look like:

```
data/raw/
├── MUTAG/MUTAG/
│   ├── MUTAG_A.txt
│   ├── MUTAG_graph_indicator.txt
│   ├── MUTAG_graph_labels.txt
│   └── MUTAG_node_labels.txt
├── NCI1/NCI1/   (same structure)
├── NCI109/NCI109/
└── ZINC/
    └── zinc250k.smi
```

> Scripts skip files that already exist. To force a fresh download, delete the target folder and re-run.

### ZINC20 full dataset (tranches)

The 250k subset above is enough for quick experiments. To reproduce the full thesis results (~16 M molecules), download ZINC20 by tranches using the URI list already included in the repo:

```console
python scripts/download/download_data/download_zinc_from_uri.py --uri scripts/download/download_data/ZINC-downloader-2D-smi.uri --workers 8
```

Downloads land in `data/raw/ZINC20/` as individual `.smi` tranche files. The script retries failed files automatically (5 attempts, exponential backoff).

**Useful options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--uri` | — | Path to the `.uri` file (required — see above) |
| `--out` | `data/raw/ZINC20/` | Output directory |
| `--workers` | `3` | Parallel download threads |
| `--n` | all | Limit to first N files (useful for testing) |
| `--retries` | `5` | Max attempts per file |

**Test with a small batch first:**

```console
python scripts/download/download_data/download_zinc_from_uri.py --uri scripts/download/download_data/ZINC-downloader-2D-smi.uri --n 10 --workers 4
```

> The `.uri` file was generated from the [ZINC20 download interface](https://zinc20.docking.org/tranches/home/) with filter: 2D → SMILES format. To update the tranche selection, download a new `.uri` file from ZINC20 and replace `scripts/download/download_data/ZINC-downloader-2D-smi.uri`.

---

## 4. Parse datasets to SMILES

The pipeline reads `.smi` files (space-separated, two columns with a header). ZINC/ZINC20 files use the header `smiles zinc_id`; benchmark datasets (MUTAG, NCI1, NCI109) use `smiles id`. The pipeline detects either column name automatically and stores the identifier as `molecule_id` in the database. MUTAG, NCI1, and NCI109 are in TU Dortmund graph format and need to be converted first; ZINC files are already in the correct format.

```console
python scripts/download/parsing_data/parse_mutag_to_smiles.py
python scripts/download/parsing_data/parse_nci1_to_smiles.py
python scripts/download/parsing_data/parse_nci109_to_smiles.py
```

This creates:

```
data/processed/
├── MUTAG/mutag_smiles.smi
├── NCI1/nci1_smiles.smi
└── NCI109/nci109_smiles.smi
```

> Parsing requires rdkit. Rows that rdkit cannot reconstruct into a valid molecule are silently skipped.

---

## 5. Run the pipeline

```console
# MUTAG
python examples/run_pipeline.py --data data/processed/MUTAG/mutag_smiles.smi

# NCI1
python examples/run_pipeline.py --data data/processed/NCI1/nci1_smiles.smi

# NCI109
python examples/run_pipeline.py --data data/processed/NCI109/nci109_smiles.smi

# ZINC 250k
python examples/run_pipeline.py --data data/raw/ZINC/zinc250k.smi

# Multiple workers for large datasets
python examples/run_pipeline.py --data data/raw/ZINC/zinc250k.smi --workers 4
```

### ZINC20 tranches (sequential run over all files)

After downloading tranches to `data/raw/ZINC20/` (see section 3), run the pipeline on every `.smi` file with:

```console
python examples/run_pipeline.py --data-dir data/raw/ZINC20 --workers 4
```

Each tranche produces its own SQLite database in `results/` (e.g. `results/AAAA.db`) and appends one summary row to the shared `results/summary.csv`.

**Run on a single `.smi` file**

To process one specific tranche, run:

```console
python examples/run_pipeline.py --data data/raw/ZINC20/AAAA.smi --workers 4
```

**If a run is interrupted**, simply re-run the command — files whose results already exist in the database are not reprocessed (the DB row for that tranche is updated in place).

---

**All CLI options:**

| Argument | Default | Description |
|---|---|---|
| `--data` | — | Path to a single `.smi` input file |
| `--data-dir` | — | Directory of `.smi` files; processes all of them in sorted order |
| `--indices` | — | Comma-separated 0-based row indices to process (e.g. `0,5,10`); only valid with `--data` |
| `--sample` | `300` | Molecules sampled to auto-estimate WL step count K |
| `--cap` | `50` | Max WL iterations allowed during K estimation |
| `--db` | `results/<dataset_stem>.db` | SQLite database path; derived from input filename by default (not valid with `--data-dir`) |
| `--summary-csv` | `results/summary.csv` | CSV with one summary row per dataset |
| `--workers` | `1` | Parallel worker processes (`multiprocessing.Pool`) |

Results are written to `results/` (created automatically, not tracked by git).

---

### Run selected molecules from an existing `.smi` file

Use `--indices` to process only specific molecules by their 0-based row index (header row not counted).

```console
python examples/run_pipeline.py --data data/processed/MUTAG/mutag_smiles.smi --indices 0,5,10 --db results/selected_molecules.db --workers 1
```

The original `.smi` file is not modified.  
Only a temporary file is created, processed, and removed after the run.

The results are saved to:

```text
results/selected_molecules.db
```
---

## 6. What the pipeline does

For each molecule:

1. **Parse SMILES** with RDKit → NetworkX graph with `atomic_num`, `is_aromatic`, `formal_charge` on nodes and `bond_type` on edges.
2. **Run topological 1-WL** for K steps (no node labels) → count color classes `n_colors_top`.
3. **Run atom-aware 1-WL** for K steps (initialized from `atomic_num`) → count color classes `n_colors_atom`.
4. **Build flip graph** from the stable atom-aware WL partition using the majority-edge rule.
5. **Check skeleton conditions** (Lemma 16, Kiefer 2020) on the flip graph.
6. **Bouquet forest check** — each connected component of the flip graph is classified as:
   - *tree* (skipped),
   - *bouquet* (C₅ with identical tree petals at each cycle vertex),
   - or *invalid* (not a bouquet).
   A molecule is marked as **non-identifiable** if every component is a tree or bouquet AND at least two bouquet components have the same petal signature.

K is estimated automatically from a random sample of `--sample` molecules before the main run.

---

## 7. Output

```
results/
├── mutag_smiles.db   # SQLite — one DB per dataset, one row per molecule
├── nci1_smiles.db
├── nci109_smiles.db
├── <tranche_name>.db # one file per ZINC20 tranche
└── summary.csv       # one row per dataset; re-running replaces the existing row
```

Key columns in the per-dataset `.db` files (e.g. `results/mutag_smiles.db`):

| Column | Description |
|--------|-------------|
| `molecule_id` | Identifier from the `.smi` file |
| `dataset_name` | Source file name (e.g. `mutag_smiles.smi`) |
| `smiles` | Input SMILES |
| `n_nodes` / `n_edges` | Graph size |
| `wl_stats` | JSON blob with per-molecule WL metrics: `n_colors_top`, `n_colors_atom`, `wl_converged_top/atom`, `wl_iters_top/atom`, `color_ratio_atom_to_top` |
| `top_bouquet_forest_verdict` | `1` = molecule is identifiable under topological WL; `0` = not identifiable or invalid flip structure |
| `atom_bouquet_forest_verdict` | `1` = molecule is identifiable under atom-aware WL; `0` = not identifiable or invalid flip structure |

**`summary.csv` columns:**

| Column | Description |
|--------|-------------|
| `dataset_name` | Source file name |
| `total_molecules` | Total molecules processed |
| `parsed_ok` / `parsed_ok_pct` | Successfully parsed count and percentage |
| `bf_verdict_top` | Identifiable count (topological WL) |
| `bf_verdict_atom` | Identifiable count (atom-aware WL) |
| `k_max` / `k_p95` | WL step estimates (max and 95th percentile) |
| `avg_colors_top` / `avg_colors_atom` | Average WL color class counts |

**Example SQLite queries:**

```sql
-- Count identifiable molecules per dataset
SELECT dataset_name,
       SUM(atom_bouquet_forest_verdict) AS identifiable,
       COUNT(*) AS total
FROM results
GROUP BY dataset_name;

-- Fetch all non-identifiable molecules from MUTAG
SELECT molecule_id, smiles, n_nodes, n_edges
FROM results
WHERE dataset_name = 'mutag_smiles.smi'
  AND atom_bouquet_forest_verdict = 0;
```

---

## 8. Reproducibility

The full results reported in the thesis (16,085,612 molecules) were produced with the following command sequence:

```console
pixi install
pixi run test
python examples/run_reproducibility.py --workers 8
```

`run_reproducibility.py` processes MUTAG, NCI1, NCI109, ZINC 250k, and all ZINC20 tranches in order. It skips files that have not been downloaded yet and prints a warning for each missing file.

**Approximate runtimes** (8-core CPU, ~16 GB RAM):

| Dataset | Molecules | Wall time |
|---------|-----------|-----------|
| MUTAG | 188 | < 1 min |
| NCI1 | 4 110 | ~2 min |
| NCI109 | 4 127 | ~2 min |
| ZINC20 | ~16 M | ~12 h |

Results are written to `results/<dataset_stem>.db` (one SQLite file per dataset) and `results/summary.csv` as the pipeline runs (crash-safe via SQLite transactions).

---

## 9. Package structure

```
src/wl_identifiability/      # Core package — import from here
├── wl.py                    # 1-WL algorithm
├── graph_construction.py    # RDKit Mol → NetworkX Graph
├── flip_graph.py            # Flip graph construction
├── skeleton.py              # Skeleton graph, Lemma 16 checks
├── bouquet.py               # Bouquet detection (Theorem 17)
├── experiments.py           # is_smi/mol/graph_identifiable, analyze_single_molecule, run_molecule_analysis_pipeline
└── visualization.py         # WL coloring visualisation

examples/
├── run_pipeline.py          # Main CLI entry point
└── run_reproducibility.py  # Cross-platform reproducibility runner

scripts/
├── database/db.py           # SQLite schema and insert functions
├── examples/run_pipeline.py # Alternative CLI with more options (--limit, --wl-steps)
└── download/                # Download and parsing scripts
```

---

## 10. Tests and developer checks

```console
# Run the test suite
pixi run test
# or: pytest tests/ -v --tb=short

# Lint (src/wl_identifiability/, tests/, scripts/, examples/)
pixi run lint

# Format (same scope as lint)
pixi run fmt

# Type-check (src/wl_identifiability/ only)
pixi run typecheck
```

The test suite covers WL coloring, graph construction, flip graph, skeleton, bouquet detection, experiments, and visualisation.

---

## 11. Troubleshooting

**`ImportError: No module named 'wl_identifiability'`**
The package is not installed. Run `pixi install` or `pip install -e .` from the project root.

**`ImportError: No module named 'rdkit'`**
rdkit is not available via pip on all platforms. Use `pixi install` or:
```console
conda install -c conda-forge rdkit
```

**`FileNotFoundError: data/processed/MUTAG/mutag_smiles.smi`**
`mutag_smiles.smi` is included in the repository — no download or parsing needed. Check that the repository was cloned correctly and that `data/processed/MUTAG/mutag_smiles.smi` exists.

**`FileNotFoundError: data/raw/MUTAG/MUTAG/MUTAG_A.txt`**
Run the download script first (step 3).

**Pipeline is slow on large datasets**
Use `--workers N` where N matches your CPU count. For ZINC 250k, `--workers 4` is a good starting point.

**SQLite DB already exists**
`INSERT OR REPLACE` is used, so re-running the same dataset overwrites existing rows. To start fresh, delete the per-dataset file (e.g. `results/mutag_smiles.db`).

---

## Citation

If you use this code or results in your work, please cite:

```bibtex
@thesis{mykytyshyn2026wl,
  author  = {Yaryna Mykytyshyn},
  title   = {Is 1-WL Expressivity Sufficient for Molecular Graphs?},
  year    = {2026},
}
```

---

## Author

**Yaryna Mykytyshyn**, 2026.
