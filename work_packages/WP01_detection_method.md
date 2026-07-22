# WP01: Improve Misfit Detection Method

## Goal

Produce a defensible diagnostic taxonomy and v1 default test set for nonlinear SEM misfit detection. This work package decides what the paper and package should test, how each test should be interpreted, and which diagnostics are too unstable or unclear for the default workflow.

## Primary Files

- `SEMmisfit.R`
- `misfit.tex`
- `SEMmisfit.bib`
- Optional new notes file: `detection_method_notes.md`

## Read First

- In `misfit.tex`, read sections `Residuals in SEM`, `Detecting Nonlinear Misfit in SEM`, and `Simulation Study - Misfit Detection`.
- In `SEMmisfit.R`, read the residual extraction block, tests applied to `stdresiduals`, the HFI block, and result summarization.
- In `SEMmisfit.bib`, inspect existing entries for HFI, dHSIC, energy statistics, and modern nonparametric independence tests.

## Do Not Change

- Do not rewrite the manuscript prose except for small notes or TODO markers if necessary.
- Do not regenerate `simtable.tex`, `sim.pdf`, or `simresults.rda`.
- Do not build the R package API here; WP02 owns that.

## Required Output

Create or update `detection_method_notes.md` with:

- A taxonomy of diagnostic targets:
  - covariance/mean misfit,
  - marginal residual distribution,
  - joint residual distribution,
  - residual dependence/nonlinear association,
  - heteroscedasticity/subgroup structure,
  - localization by variable, pair, case, or moderator.
- A comparison of residual representations:
  - whitened observed residuals from implied mean/covariance,
  - conditional residuals for each variable given the rest,
  - rowwise squared residual distances,
  - pairwise residuals.
- A calibration plan:
  - naive asymptotic p-values,
  - permutation tests,
  - parametric bootstrap from the fitted SEM,
  - sample splitting or cross-fitting if useful.
- A recommended v1 default diagnostic set, with each test mapped to its null hypothesis and interpretation.
- A short list of diagnostics to exclude from defaults unless later evidence supports them.

## Acceptance Criteria

- Every recommended diagnostic has a clear null hypothesis.
- The notes explicitly distinguish distributional diagnostics from residual-dependence diagnostics.
- The notes state how fitted-parameter uncertainty affects p-value interpretation.
- The notes identify what WP02 should expose as package defaults.
- The notes identify what WP05 must simulate.

## Handoff Notes

After finishing, update `CURRENT_HANDOFF.md` with:

- the chosen default diagnostics,
- unresolved calibration concerns,
- any manuscript claims that must be softened,
- exact next recommended work package.

