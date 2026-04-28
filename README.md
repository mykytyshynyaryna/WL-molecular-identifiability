# Is 1-WL Expressivity Sufficient for Molecular Graphs?

This project analyses whether the **1-Weisfeiler-Lehman (1-WL) graph isomorphism test** can uniquely identify every molecule in standard benchmark datasets. The pipeline runs WL coloring on molecular graphs, builds a flip graph from the stable color partition, and checks each connected component against the bouquet forest criterion (Theorem 17, Kiefer 2020) to detect non-identifiable molecules.

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

**Alternative (pip):**

```bash
pip install -e ".[dev]"
```

> If rdkit fails on pip, install it via conda first:
> ```bash
> conda install -c conda-forge rdkit
> pip install -e ".[dev]" --no-deps
> ```

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
│   ├── MUTAG_A.txt              # edge list
│   ├── MUTAG_graph_indicator.txt
│   ├── MUTAG_graph_labels.txt
│   └── MUTAG_node_labels.txt
├── NCI1/NCI1/   (same structure)
├── NCI109/NCI109/
└── ZINC/
    └── zinc250k.smi             # already in SMILES format
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

Each `.smi` file has two columns with a header:

```
smiles zinc_id
CC1=CC2=C(...)  1
O=C(...)        2
```

The `zinc_id` column holds the original graph index (1-based) from the TU dataset.

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

Terminal output at the end of each run:

```
=== Pipeline summary ===
Dataset           : data/processed/MUTAG/mutag_smiles.smi
Sample size       : 300
...
Total molecules   :      188
Parsed OK         :      188  (100.0%)
Parse/WL errors   :        0
BF verdict (atom) :      188
```

---

## 7. Package structure

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
└── download/                # Download and parsing scripts
```

---

## 8. Troubleshooting

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
