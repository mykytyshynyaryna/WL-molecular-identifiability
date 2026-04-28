
### Summary of bouquet-forest verdicts (16,085,612 molecules across MUTAG, NCI1, NCI109, ZINC20)

- **`top=1, atom=0` (suspicious):** 0 (0.000000%)
- **`top=0, atom=1` (expected mismatch):** 62,988 (0.3916%)
- **`top=1, atom=1` (both accept):** 15,990,031 (99.41%)
- **`top=0, atom=0` (both reject):** 32,591 (0.2026%)

**Top-level rejection reasons (raw):**
- `cycle_not_5`: 58,242
- `petals_not_isomorphic`: 34,089
- `multiple_cycles_in_component`: 3,248

**Atom-level rejection reasons (raw):**
- `cycle_not_5`: 28,997
- `petals_not_isomorphic`: 2,416
- `multiple_cycles_in_component`: 1,178

**Reason categories (mapped):**

```
                         top   atom
category                           
not_c5_cycle           58242  28997
petals_not_isomorphic  34089   2416
multiple_cycles         3248   1178
```

**Representative-sample diversity (by `dataset_name`):**

- `top1_atom0`: 0 rows
- `top0_atom1`: 12 rows across 6 datasets (ADAD.smi, BBAD.smi, BCAD.smi, BDAD.smi, BEAD.smi, CBAD.smi)
- `top1_atom1`: 12 rows across 6 datasets (BBAD.smi, BCAD.smi, BDAD.smi, BDED.smi, BEAD.smi, CBAD.smi)
- `top0_atom0`: 12 rows across 6 datasets (ADAD.smi, BCAD.smi, BDAD.smi, BDED.smi, BEAD.smi, CBAD.smi)

**Qualitative observations:**

- The `top1_atom0` set is the critical pattern for the thesis. Empirically it is
  **empty** in this DB —
  no molecule passes the coarse top-level test while failing the finer atom-level one.
- `top0_atom1` (n=62,988) is the inverse: top-level rejects, atom-level accepts.
  Structurally these are cases where the coarse view is more conservative (e.g. petals
  collapse at the atom level into structures that satisfy the bouquet definition).
- The dominant rejection cause is `cycle_not_5` — molecules whose cyclic component is not
  a 5-cycle. The next is `petals_not_isomorphic` (the cycle exists but the dangling
  attached subgraphs are not all isomorphic), and the rarest is
  `multiple_cycles_in_component`.

Outputs saved to `C:\Users\yaryn\OneDrive\Desktop\Thesis\Is_1_WL_Expressivity_Sufficient_for_Molecular_Graphs\scripts\results\interesting_cases`:
- `verdict_pattern_counts_by_dataset.csv`
- `sample_<pattern>.csv` (×4) and `grid_<pattern>.svg` (+ `.pdf` when `cairosvg` is installed) (where the subset is non-empty)
- `sample_reason_<category>.csv` and `grid_reason_<category>.svg` (+ `.pdf` when `cairosvg` is installed) per non-empty category
- `top_reasons_raw.csv`, `atom_reasons_raw.csv`, `reason_categories.csv`
- `verdict_disagreements.csv`
- `summary.md` (this block)
