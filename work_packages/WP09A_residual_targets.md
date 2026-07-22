# WP09A: Residual Targets And Global Representations

## Goal

Define the diagnostic nulls and residual representations used by the revised
prototype. Replace Cholesky whitening as the default basis for global checks
while retaining it as a clearly labelled optional comparison.

## Primary Files

- `sem_misfit_api_prototype.R`
- `detection_method_notes.md`
- `misfit.tex`
- `CURRENT_HANDOFF.md`

## Read First

- `CURRENT_HANDOFF.md`
- `detection_method_notes.md`
- In `sem_misfit_api_prototype.R`, read `compute_sem_residuals()`,
  `compute_sem_diagnostics()`, and every use of `whitened`.
- In `misfit.tex`, read the residual-representation and prototype-software
  subsections.

## Do Not Change

- Do not implement conditional-residual tests or plots beyond preserving the
  existing residual calculation; WP09B owns localization.
- Do not implement bootstrap calibration or regenerate simulation outputs.
- Do not redesign the empirical example or package structure.

## Required Output

- State exactly three targets in code labels, notes, and manuscript prose:
  mean/covariance fit; full Gaussian observed-data fit; and conditional
  adequacy for a named variable.
- Implement numerically safeguarded symmetric/ZCA whitening
  $Z=(X-\hat\mu)\hat\Sigma^{-1/2}$ as the default global representation.
- Preserve the existing Cholesky transform as `cholesky_innovations`, with an
  explicit order-dependence warning and no original-variable pair labels.
- Ensure global distribution, global dHSIC, marginal distribution, and
  case-distance diagnostics state which representation they use.
- Revise notes and manuscript wording so a rejection reports its exact null,
  rather than generic nonlinear SEM misfit.

## Acceptance Criteria

- Symmetric-whitened residual covariance is numerically identity for a
  positive-definite non-diagonal implied covariance.
- Reordering input variables reorders symmetric-whitened output without
  changing its values apart from the same permutation.
- Cholesky coordinate order sensitivity is demonstrated and documented.
- No pairwise output or manuscript text equates a whitened coordinate with an
  original-variable residual.

## Handoff Notes

Record the selected eigenvalue tolerance, affected API labels, manuscript
sections changed, test evidence, and whether WP09B may begin.
