# WP09B: Conditional-Residual Localization

## Goal

Make variable-specific conditional residuals the default localization object
and test their behavior independently of global whitening diagnostics.

## Primary Files

- `sem_misfit_api_prototype.R`
- New or existing prototype test file(s)
- `detection_method_notes.md`
- `misfit.tex`

## Read First

- `CURRENT_HANDOFF.md`
- `work_packages/WP09A_residual_targets.md`
- The completed WP09A handoff entry.
- Existing `compute_conditional_residuals()`, pairwise output, and conditional
  plot code.

## Do Not Change

- Do not alter symmetric-whitening definitions settled in WP09A.
- Do not implement refit bootstrap calibration or regenerate simulations.
- Do not make conditional residuals mutually independent by construction or
  test them against each other as a global independence family.

## Required Output

- Implement variable-level tests of $R_j \perp X_{-j}$ and directional pair
  tests of $R_j \perp X_k$ using variable-labelled conditional residuals.
- Add optional conditional-scale diagnostics of $R_j^2$ against $X_{-j}$.
- Return separate `conditional_variables`, `conditional_pairs`, and
  `conditional_scale` output tables with statistic, p-value/effect size,
  calibration label, and exact null.
- Update conditional-residual plots to support predictors, fitted values, and
  moderators using original variable labels.
- Remove symmetric- or Cholesky-whitened pair tables from default substantive
  localization. They may remain as explicitly transformed exploratory output.

## Acceptance Criteria

- Conditional residuals agree with direct Gaussian linear conditional
  regressions on simulated non-diagonal covariance data.
- Directional output distinguishes `R_j` versus `X_k` from `R_k` versus
  `X_j`.
- Tests preserve variable labels and handle near-singular conditioning blocks
  with informative errors or documented regularization.
- The manuscript states that a conditional rejection can reflect mean, scale,
  or broader conditional-distribution inadequacy, not uniquely a nonlinear
  mean term.

## Handoff Notes

Record the chosen conditional tests, multiplicity strategy placeholder,
numerical handling, tests run, and the interfaces WP10A must calibrate.
