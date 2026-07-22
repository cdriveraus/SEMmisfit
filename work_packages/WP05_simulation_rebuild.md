# WP05: Simulation Rebuild And Output Alignment

## Goal

Align `SEMmisfit.R`, regenerated outputs, and manuscript claims after WP01 settles the diagnostic set. Simulations should illustrate diagnostic separation rather than claim a definitive benchmark.

## Primary Files

- `SEMmisfit.R`
- `simresults.rda`
- `simtable.tex`
- `sim.pdf`
- `misfit.tex`

## Read First

- `CURRENT_HANDOFF.md`
- `detection_method_notes.md`
- In `SEMmisfit.R`, read simulation conditions, model specification, residual extraction, diagnostic tests, summarization, table output, and figure output.
- In `misfit.tex`, read `Simulation Study - Misfit Detection` and `Simulation Results`.

## Do Not Change

- Do not change the manuscript framing outside the simulation section unless needed for consistency.
- Do not add new diagnostics not approved by WP01.
- Do not treat exploratory results as final evidence without checking calibration.

## Required Output

Revise and rerun the simulation workflow so that:

- Conditions map cleanly to diagnostic targets:
  - no misfit,
  - omitted covariance/linear relation,
  - omitted interaction,
  - nonlinear monotone or nonmonotone relation if useful,
  - skewed marginal distribution without residual dependence,
  - heavy tails or outliers,
  - group mean/variance heterogeneity.
- Diagnostics match WP01 defaults.
- Output table and figure use clear labels and stable formatting.
- Manuscript text describes the simulation design and results accurately.

## Acceptance Criteria

- `SEMmisfit.R` can regenerate `simtable.tex` and `sim.pdf`.
- Detection rates in the manuscript match generated outputs.
- False-positive behavior under no-misfit and negative-control conditions is discussed.
- Claims are phrased as diagnostic illustrations, not universal performance rankings.
- Any long runtime or package dependency issue is documented in `CURRENT_HANDOFF.md`.

## Handoff Notes

Update `CURRENT_HANDOFF.md` with:

- command used to rerun simulations,
- runtime and package issues,
- generated files,
- key result patterns,
- remaining simulation limitations.

