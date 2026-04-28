# Dataset Download Scripts

All scripts live in `scripts/download/`.  
Run them from the **project root** directory.

---

## Quick start — download everything

```bash
python scripts/download/download_all.py
```

Downloads all six datasets in sequence and prints a summary at the end.

---

## Individual dataset scripts

| Script | Dataset | Source | Format |
|--------|---------|--------|--------|
| `download_mutag.py` | MUTAG | TU Dortmund | TU plain-text (ZIP) |
| `download_nci1.py` | NCI1 | TU Dortmund | TU plain-text (ZIP) |
| `download_nci109.py` | NCI109 | TU Dortmund | TU plain-text (ZIP) |
| `download_dd.py` | D&D (DD) | TU Dortmund | TU plain-text (ZIP) |
| `download_enzymes.py` | ENZYMES | TU Dortmund | TU plain-text (ZIP) |
| `download_zinc.py` | ZINC 250k | chemical_vae / GitHub | `.smi` |
| `download_zinc_from_uri.py` | ZINC (custom tranche) | ZINC15 URI list | `.smi` |

Run any single script directly:

```bash
python scripts/download/download_mutag.py
python scripts/download/download_nci1.py
python scripts/download/download_nci109.py
python scripts/download/download_dd.py
python scripts/download/download_enzymes.py
python scripts/download/download_zinc.py

# Tranche downloader (requires a local .uri file from ZINC15)
python scripts/download/download_zinc_from_uri.py --uri path/to/file.uri
```

---

## Output layout

After running, `data/raw/` will contain:

```
data/raw/
├── MUTAG/
│   ├── MUTAG.zip
│   └── MUTAG/
│       ├── MUTAG_A.txt
│       ├── MUTAG_graph_indicator.txt
│       ├── MUTAG_graph_labels.txt
│       ├── MUTAG_node_labels.txt
│       └── README.txt
├── NCI1/
│   ├── NCI1.zip
│   └── NCI1/  (same structure)
├── NCI109/
│   ├── NCI109.zip
│   └── NCI109/
├── DD/
│   ├── DD.zip
│   └── DD/
├── ENZYMES/
│   ├── ENZYMES.zip
│   └── ENZYMES/
│       └── ...  (includes ENZYMES_node_attributes.txt)
└── ZINC/
    └── zinc250k.smi   ← one SMILES per line
```

---

## Dataset notes

### TU Dortmund format (MUTAG / NCI1 / NCI109 / DD / ENZYMES)

The TU Dortmund benchmark collection distributes graph datasets as plain-text
files.  Key files:

| File | Content |
|------|---------|
| `<NAME>_A.txt` | Edge list: one `u, v` pair per line (1-indexed) |
| `<NAME>_graph_indicator.txt` | Graph ID for each node (1-indexed) |
| `<NAME>_graph_labels.txt` | Binary or multi-class label per graph |
| `<NAME>_node_labels.txt` | Integer label per node (atom/residue type) |
| `<NAME>_node_attributes.txt` | Continuous features per node (ENZYMES only) |

> **Note on SDF format:**  NCI1 and NCI109 were originally distributed by the
> NCI DTP program in `.sdf` format.  If you need SDF files, see the `# TODO`
> comments inside `download_nci1.py` and `download_nci109.py` for the NCI DTP
> download portal URL.

### ZINC (zinc250k.smi)

The ZINC250k subset contains 250 000 drug-like SMILES strings drawn from the
ZINC15 database.  The file has one SMILES string per line with no header.

This is the same subset used in the chemical VAE (Gómez-Bombarelli et al. 2018)
and the junction tree VAE (Jin et al. 2018) papers.

If you need a larger or more specific ZINC subset, use the ZINC15 tranche
download interface at <https://zinc15.docking.org/tranches/home/> and select
"2D" → "SMILES" format.  Save the resulting URI list and run
`download_zinc_from_uri.py` with the `--uri` argument.

---

## Requirements

All scripts use the **Python standard library only** (`urllib`, `zipfile`,
`csv`, `pathlib`, `subprocess`).  No third-party packages are needed.

Python ≥ 3.8 required.

---

## Re-downloading

Scripts skip files that already exist.  To force a fresh download, delete the
target directory (e.g. `data/raw/MUTAG/`) or the specific file, then re-run.
