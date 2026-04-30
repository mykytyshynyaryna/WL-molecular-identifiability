# Is 1-WL Expressivity Sufficient for Molecular Graphs?

> **Two branches in this repository:**
> - [`main`](../../tree/main) — implementation only: the core `wl_identifiability` package, pipeline scripts, tests, and download utilities. No analysis, no results.
> - [`thesis`](../../tree/thesis) *(this branch)* — everything in `main`, plus: research results across 16 M molecules, Jupyter notebooks with per-dataset analysis, profiling scripts, and the full thesis writeup context.

This project investigates whether the **1-Weisfeiler-Lehman (1-WL) graph isomorphism test** can uniquely identify every molecule in standard benchmark datasets. The pipeline runs WL coloring on molecular graphs, builds a flip graph from the stable color partition, and checks each connected component against the bouquet forest criterion to detect non-identifiable molecules.

---

## Background

**1-Weisfeiler-Lehman (1-WL) test** is an iterative colour-refinement algorithm on graphs. Each node starts with an initial colour (label); at every step every node aggregates the multiset of its neighbours' colours and hashes it into a new colour. The process stops when no colour changes. Two graphs are *non-isomorphic* by 1-WL if they produce different final colour histograms, and *potentially isomorphic* otherwise.

**Two WL modes used in this project:**
- *Topological* — initial node colour is uniform (structure only, no atom labels).
- *Atom-aware* — initial node colour encodes the atomic symbol, so chemical identity is used from the first iteration.

**Bouquet forest** (Kiefer 2020, Definition 9 / Theorem 17): a graph is a *bouquet* if it consists of a single C₅ cycle with zero or more pendant trees ("petals") attached at one vertex. A *bouquet forest* is a disjoint union of isomorphic bouquets. Theorem 17 states that a graph is *not* uniquely identified by 1-WL if and only if its flip graph decomposes into a bouquet forest with at least two components.

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

**Interpreting verdicts:** A molecule is called *identifiable* (verdict = 1) when its flip graph is a bouquet forest — meaning 1-WL is sufficient to distinguish it from all structurally similar molecules. "Both WL modes accept" means identifiable under both topological and atom-aware coloring. "Both reject" means 1-WL cannot distinguish the molecule regardless of whether atom labels are used.

---

## Prerequisites

- Python 3.10–3.12
- [pixi](https://pixi.sh) (recommended — handles rdkit automatically via conda-forge)

---

## 1. Install

**Recommended:**

```bash
pixi install
```

This creates an isolated environment with all dependencies (networkx, numpy, polars, rdkit, matplotlib, pytest).

**Key dependency versions (resolved by pixi):**

| Package | Version |
|---------|---------|
| Python | 3.10 – 3.12 |
| RDKit | ≥ 2023.09 |
| NetworkX | ≥ 3.1 |
| pandas | ≥ 2.0 |
| pytest | ≥ 7.0 |

**Alternative (pip):**

```bash
pip install -e ".[dev]"
```

> If rdkit fails on pip, install it via conda first:
> ```bash
> conda install -c conda-forge rdkit
> pip install -e ".[dev]" --no-deps
> ```

**Recommended first steps after install:**

```bash
# 1. Verify installation
pixi run test

# 2. Quick smoke test
python examples/run_pipeline.py --data data/processed/MUTAG/mutag_smiles.smi --sample 50
```

---

## 2. Download datasets

All datasets go into `data/raw/`. Run from the project root.

**Download everything at once:**

```bash
python "scripts/download/download data/download_all.py"
```

**Or download individually:**

```bash
python "scripts/download/download data/download_mutag.py"
python "scripts/download/download data/download_nci1.py"
python "scripts/download/download data/download_nci109.py"
python "scripts/download/download data/download_zinc.py"
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

---

## 3. Parse datasets to SMILES

The pipeline reads `.smi` files (space-separated SMILES with an ID column). ZINC is already in this format. MUTAG, NCI1, and NCI109 are in TU Dortmund graph format and need to be converted first.

```bash
python "scripts/download/parsing data/parse_mutag_to_smiles.py"
python "scripts/download/parsing data/parse_nci1_to_smiles.py"
python "scripts/download/parsing data/parse_nci109_to_smiles.py"
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

## 4. Run the pipeline

```bash
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

**All CLI options:**

| Argument | Default | Description |
|---|---|---|
| `--data` | `data/AAAA.smi` | Path to `.smi` input file |
| `--sample` | `300` | Molecules sampled to auto-estimate WL step count K |
| `--cap` | `50` | Max WL iterations allowed during K estimation |
| `--db` | `results/all_results.db` | SQLite database — all datasets share one file |
| `--summary-csv` | `results/summary.csv` | CSV with one summary row per dataset |
| `--workers` | `1` | Parallel worker processes (`multiprocessing.Pool`) |

Results are written to `results/` (created automatically, not tracked by git).

---

## 5. What the pipeline does

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

## 6. Output

```
results/
├── all_results.db    # SQLite — one row per molecule, tagged by dataset_name
└── summary.csv       # one row per dataset; re-running replaces the existing row
```

Key columns in `all_results.db`:

| Column | Description |
|--------|-------------|
| `molecule_id` | Identifier from the `.smi` file |
| `dataset_name` | Source file name (e.g. `mutag_smiles.smi`) |
| `smiles` | Input SMILES |
| `n_nodes` / `n_edges` | Graph size |
| `n_colors_top` / `n_colors_atom` | WL color counts (topological / atom-aware) |
| `top_bouquet_forest_verdict` | `1` = flip graph is a bouquet forest (topological WL) |
| `atom_bouquet_forest_verdict` | `1` = flip graph is a bouquet forest (atom-aware WL) |

**`summary.csv` columns:**

| Column | Description |
|--------|-------------|
| `dataset_name` | Source file name |
| `total` | Total molecules processed |
| `identifiable_top` | Identifiable count (topological WL) |
| `identifiable_atom` | Identifiable count (atom-aware WL) |
| `pct_identifiable_atom` | Percentage identifiable (atom-aware) |

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

## 7. Notebooks

Jupyter notebooks in `notebooks/` contain per-dataset analysis of non-identifiable molecules.

| Notebook | Contents |
|----------|----------|
| `interesting_cases_analysis_MUTAG.ipynb` | Size/degree distribution of non-identifiable MUTAG molecules; visualisation of WL coloring on representative examples |
| `interesting_cases_analysis_NCI1.ipynb` | Breakdown of rejection reasons for NCI1; structural motifs common in non-identifiable anti-cancer compounds |
| `interesting_cases_analysis_NCI109.ipynb` | Comparison of NCI109 rejection patterns vs NCI1; overlap between the two datasets |
| `interesting_cases_analysis_ZINC20.ipynb` | Scale analysis of 16 M ZINC20 molecules; histogram of non-identifiable SMILES and their pharmacophore classes |

```bash
jupyter notebook notebooks/
```

---

## 8. Reproducibility

The full results reported in the thesis (16,085,612 molecules) were produced with the following command sequence:

```bash
pixi install
pixi run test

for f in data/processed/*.smi data/raw/ZINC/zinc250k.smi; do
    python examples/run_pipeline.py --data "$f" --workers 8
done
```

**Approximate runtimes** (8-core CPU, ~16 GB RAM):

| Dataset | Molecules | Wall time |
|---------|-----------|-----------|
| MUTAG | 188 | < 1 min |
| NCI1 | 4 110 | ~2 min |
| NCI109 | 4 127 | ~2 min |
| ZINC20 | ~16 M | ~12 h |

Results are written to `results/all_results.db` and `results/summary.csv` as the pipeline runs (crash-safe via SQLite transactions).

---

## 9. Package structure

```
src/wl_identifiability/      # Core package — import from here
├── wl.py                    # 1-WL algorithm
├── graph_construction.py    # RDKit Mol → NetworkX Graph
├── flip_graph.py            # Flip graph construction
├── skeleton.py              # Skeleton graph, Lemma 16 checks
├── bouquet.py               # Bouquet detection (Theorem 17)
├── experiments.py           # analyze_single_molecule, run_molecule_analysis_pipeline
└── visualization.py         # WL coloring visualisation

examples/
└── run_pipeline.py          # Main CLI entry point

scripts/
├── database/db.py           # SQLite schema and insert functions
├── pipeline/run_pipeline.py # Alternative CLI with more options (--limit, --wl-steps)
├── benchmarks/              # Flip graph benchmarks
├── profiling/               # cProfile scripts
└── download/                # Download and parsing scripts
```

---

## 10. Tests

```bash
pixi run test
# or
pytest tests/ -v --tb=short
```

The suite covers WL coloring, graph construction, flip graph, skeleton, bouquet detection, experiments, and visualisation.

---

## 11. Troubleshooting

**`ImportError: No module named 'wl_identifiability'`**
The package is not installed. Run `pixi install` or `pip install -e .` from the project root.

**`ImportError: No module named 'rdkit'`**
rdkit is not available via pip on all platforms. Use `pixi install` or:
```bash
conda install -c conda-forge rdkit
```

**`FileNotFoundError: data/processed/MUTAG/mutag_smiles.smi`**
Run the parsing script first (step 3).

**`FileNotFoundError: data/raw/MUTAG/MUTAG/MUTAG_A.txt`**
Run the download script first (step 2).

**Pipeline is slow on large datasets**
Use `--workers N` where N matches your CPU count. For ZINC 250k, `--workers 4` is a good starting point.

**SQLite DB already exists**
`INSERT OR REPLACE` is used, so re-running the same dataset overwrites existing rows. To start fresh, delete `results/all_results.db`.

---

## Citation

If you use this code or results in your work, please cite:

```bibtex
@mastersthesis{mykytyshyn2026wl,
  author  = {Yaryna Mykytyshyn},
  title   = {Is 1-WL Expressivity Sufficient for Molecular Graphs?},
  year    = {2026},
}
```

---

## Author

**Yaryna Mykytyshyn**, 2026.
