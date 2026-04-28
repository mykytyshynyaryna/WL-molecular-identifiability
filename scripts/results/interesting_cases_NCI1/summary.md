
### Summary of bouquet-forest verdicts (3,999 molecules across MUTAG, NCI1, NCI109, ZINC20)

- **`top=1, atom=0` (suspicious):** 0 (0.000000%)
- **`top=0, atom=1` (expected mismatch):** 54 (1.3503%)
- **`top=1, atom=1` (both accept):** 3,449 (86.25%)
- **`top=0, atom=0` (both reject):** 496 (12.4031%)

**Top-level rejection reasons (raw):**
- `cycle_not_5`: 301
- `multiple_cycles_in_component`: 224
- `petals_not_isomorphic`: 25

**Atom-level rejection reasons (raw):**
- `cycle_not_5`: 256
- `multiple_cycles_in_component`: 217
- `petals_not_isomorphic`: 23

**Reason categories (mapped):**

```
                       top  atom
category                        
not_c5_cycle           301   256
multiple_cycles        224   217
petals_not_isomorphic   25    23
```

**Representative-sample diversity (by `dataset_name`):**

- `top1_atom0`: 0 rows
- `top0_atom1`: 2 rows across 1 datasets (nci1_smiles.smi)
- `top1_atom1`: 2 rows across 1 datasets (nci1_smiles.smi)
- `top0_atom0`: 2 rows across 1 datasets (nci1_smiles.smi)

**Qualitative observations:**

- The `top1_atom0` set is the critical pattern for the thesis. Empirically it is
  **empty** in this DB —
  no molecule passes the coarse top-level test while failing the finer atom-level one.
- `top0_atom1` (n=54) is the inverse: top-level rejects, atom-level accepts.
  Structurally these are cases where the coarse view is more conservative (e.g. petals
  collapse at the atom level into structures that satisfy the bouquet definition).
- The dominant rejection cause is `cycle_not_5` — molecules whose cyclic component is not
  a 5-cycle. The next is `petals_not_isomorphic` (the cycle exists but the dangling
  attached subgraphs are not all isomorphic), and the rarest is
  `multiple_cycles_in_component`.

Outputs saved to `C:\Users\yaryn\OneDrive\Desktop\Thesis\Is_1_WL_Expressivity_Sufficient_for_Molecular_Graphs\scripts\results\interesting_cases_NCI1`:
- `verdict_pattern_counts_by_dataset.csv`
- `sample_<pattern>.csv` (×4) and `grid_<pattern>.svg` (+ `.pdf` when `cairosvg` is installed) (where the subset is non-empty)
- `sample_reason_<category>.csv` and `grid_reason_<category>.svg` (+ `.pdf` when `cairosvg` is installed) per non-empty category
- `top_reasons_raw.csv`, `atom_reasons_raw.csv`, `reason_categories.csv`
- `verdict_disagreements.csv`
- `summary.md` (this block)
