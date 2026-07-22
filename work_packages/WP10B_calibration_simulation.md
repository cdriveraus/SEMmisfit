# WP10B: Calibration Simulation And Generated Outputs

## Goal

Validate the revised diagnostic representations and calibrations with a
simulation designed for Type-I error, power, localization, and computational
cost rather than an illustrative only comparison.

## Primary Files

- `SEMmisfit.R`
- `simresults.rda`
- `simtable.tex`
- `sim.pdf`
- `misfit.tex`

## Read First

- `CURRENT_HANDOFF.md`
- `work_packages/WP10A_bootstrap_calibration.md`
- Completed WP09A, WP09B, and WP10A handoff entries.
- Current simulation conditions and result summaries in `SEMmisfit.R`.

## Do Not Change

- Do not change the residual definitions or calibration contracts chosen in
  WP09/WP10A.
- Do not rewrite the manuscript outside the simulation sections except for
  strictly necessary terminology consistency.
- Do not use a lower replication count merely to reduce runtime.

## Required Output

- Simulate Gaussian null, marginal-skewness-only, heavy-tail-only, covariance
  misspecification, nonlinear conditional mean, conditional heteroscedasticity,
  subgroup heterogeneity, rotated independent non-Gaussian disturbances, and
  higher-order-only dependence.
- Cross $n=100,500$ with 4 and 8 observed variables. Use 1,000 outer
  replications for the three null/negative-control conditions and 500 for all
  non-null conditions.
- For every eligible diagnostic, compare conditional permutation,
  parametric-fixed, and parametric-refit calibration. Use 999 replicates for
  each calibration distribution.
- Use max-statistic bootstrap calibration for variable and pairwise families;
  do not use coarse Bonferroni-adjusted permutation p-values as the primary
  localization comparison.
- Produce tables/figures with rejection rate, binomial Monte Carlo interval,
  median runtime, bootstrap failure rate, and localization accuracy where a
  true affected variable/pair exists.
- Update simulation prose with results generated from the saved artifact; call
  conservative tests conservative and limit claims to the studied conditions.

## Acceptance Criteria

- Script, saved results, tables, figures, and manuscript agree on conditions,
  replications, calibration modes, and rejection rates.
- Null behavior is reported separately for Gaussian and non-Gaussian negative
  controls.
- Any unavailable or failed calibration is clearly distinguishable from a
  non-rejection.
- A rerun command, package versions, seed strategy, runtime, and output paths
  are recorded in `CURRENT_HANDOFF.md`.

## Handoff Notes

Record completed conditions, run command, wall-clock runtime, failures,
generated files, result patterns, remaining computational limitations, and
whether WP11A may begin.
