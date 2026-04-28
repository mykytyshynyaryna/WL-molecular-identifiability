# Profiling

Compares baseline vs optimised implementations of the two main algorithmic bottlenecks in the WL identifiability pipeline.

All commands run from the **project root**.

---

## Structure

```
scripts/profiling/run all profiles/
├── profile_bouquet_baseline.py    # bouquet detection — baseline (nx.is_isomorphic)
├── profile_bouquet_optimized.py   # bouquet detection — optimised (AHU signatures)
├── profile_flip_graph_baseline.py # flip graph — dense NumPy submatrix approach
├── profile_flip_graph_optimized.py# flip graph — sparse edge-list approach
├── run_all_profiles.py            # run all four scenarios, print comparison table
├── profiling_utils.py             # shared: timing, cProfile, table rendering
└── dataset_loader.py              # shared: .smi file → (mol_id, G, labels) cases
```

Results are saved to `scripts/profiling/run all profiles/results/` (`.prof` binary + `.txt` summary per scenario).

---

## Implementations compared

### Bouquet detection

| Script | Approach |
|--------|---------|
| `profile_bouquet_baseline.py` | `nx.cycle_basis` + `nx.is_isomorphic` for pairwise petal comparison |
| `profile_bouquet_optimized.py` | Leaf-stripping C₅ finder + AHU canonical string signatures; no graph copies |

### Flip graph construction

| Script | Approach |
|--------|---------|
| `profile_flip_graph_baseline.py` | Dense `N×N` NumPy matrix, submatrix blocks via `np.ix_` |
| `profile_flip_graph_optimized.py` | Pure edge-list accumulation; no dense matrix |

---

## Run all scenarios

```bash
python "scripts/profiling/run all profiles/run_all_profiles.py"
python "scripts/profiling/run all profiles/run_all_profiles.py" --data data/AAAC.smi --max-mols 500 --reps 10
python "scripts/profiling/run all profiles/run_all_profiles.py" --no-save
```

Prints a summary table at the end:

```
=== Summary ===
+----------------------+---------+---------+------+--------+
| scenario             | total_s | ms/call | reps | n_mols |
+----------------------+---------+---------+------+--------+
| bouquet_baseline     | 3.1200  | 6.2400  | 5    | 100    |
| bouquet_optimized    | 0.4800  | 0.9600  | 5    | 100    |
| flip_graph_baseline  | 1.2300  | 2.4600  | 5    | 100    |
| flip_graph_optimized | 0.3100  | 0.6200  | 5    | 100    |
+----------------------+---------+---------+------+--------+
```

---

## Run individual scripts

```bash
python "scripts/profiling/run all profiles/profile_bouquet_baseline.py"
python "scripts/profiling/run all profiles/profile_bouquet_optimized.py"
python "scripts/profiling/run all profiles/profile_flip_graph_baseline.py"
python "scripts/profiling/run all profiles/profile_flip_graph_optimized.py"
```

### Common arguments

| Argument | Default | Description |
|---|---|---|
| `--data` | `data/AAAA.smi` | Path to `.smi` dataset file |
| `--max-mols` | all | Stop after N molecules |
| `--reps` | `5` | Timing repetitions |
| `--out` | `scripts/profiling/run all profiles/results` | Output directory |
| `--no-save` | off | Skip saving `.prof`/`.txt` files |
| `--top-n` | `20` | Slowest functions to show |

Example:

```bash
python "scripts/profiling/run all profiles/profile_bouquet_optimized.py" --data data/AAAC.smi --max-mols 500 --reps 10
```

---

## Saved outputs

| File | Description |
|------|-------------|
| `<scenario>.prof` | Binary cProfile dump — open with `snakeviz <file>.prof` |
| `<scenario>.txt` | Human-readable top-N functions report |
