# WP10A: Estimation-Aware Calibration

## Goal

Implement three explicitly distinct calibration paths for the revised
diagnostic statistics, with refit parametric bootstrap as the inferential
reference implementation.

## Primary Files

- `sem_misfit_api_prototype.R`
- New or existing prototype test file(s)
- `detection_method_notes.md`
- `CURRENT_HANDOFF.md`

## Read First

- `CURRENT_HANDOFF.md`
- `work_packages/WP09A_residual_targets.md`
- `work_packages/WP09B_conditional_localization.md`
- Completed WP09A/WP09B handoff entries.
- Existing `calibration`, permutation, and SEM-reader code.

## Do Not Change

- Do not alter the residual objects or localization targets decided in WP09.
- Do not rerun the full validation simulation; WP10B owns generated outputs.
- Do not claim OpenMx refit support unless it is actually implemented and
  tested.

## Required Output

- Expose `conditional_permutation`, `parametric_fixed`, and
  `parametric_refit` calibration modes. Preserve a clear legacy alias only if
  needed to avoid breaking current prototype calls.
- Define `parametric_fixed` as simulation from the fitted expectation with
  diagnostics recomputed using the original fitted expectation.
- Define `parametric_refit` as simulation from the fitted model, refitting the
  identical model restrictions, then recomputing diagnostics on the refit.
- Implement refitting first for lavaan and for simulation-native generic
  expectations where a refit callback is explicitly supplied.
- Calibrate raw statistics by bootstrap rank with $B=999$ by default. Do not
  run nested permutation tests inside bootstrap replicates.
- Calibrate maximum statistics across every selected variable- or pair-level
  family in each replicate for family-wise inference.
- Return successful-replicate count, failure reasons, seed, runtime, and Monte
  Carlo p-value along with the observed statistic.

## Acceptance Criteria

- A correctly specified lavaan Gaussian model can complete fixed and refit
  calibration with reproducible results under a fixed seed.
- Failed/nonconverged/improper refits are retained in a failure log and
  excluded only under an explicit documented rule.
- The reported p-value equals the bootstrap-rank calculation from retained
  statistics.
- Conditional permutation remains available and is visibly labelled as
  conditional exploratory evidence.

## Handoff Notes

Record the refit contract, supported sources, default `B`, failure rule,
runtime observations, tests run, and the exact statistics ready for WP10B.
