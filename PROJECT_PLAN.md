# SEM Misfit Project Plan

## Purpose

This project develops a manuscript and small R software package for residual-based diagnostics of nonlinear or otherwise hidden misfit in structural equation models. The goal is not to replace classical SEM fit statistics. The goal is to add a diagnostic layer that asks whether the model-implied expectations are adequate for the observed data beyond mean and covariance fit.

The manuscript should eventually argue that SEM fit assessment benefits from separating several diagnostic questions:

- Does the model reproduce the covariance and mean structure?
- Do standardized residuals match the assumed marginal or joint reference distribution?
- Are residual components mutually independent after whitening?
- Is residual misfit localized to particular variables, cases, pairs, or moderators?
- Is there evidence for heterogeneity, subgroup structure, interactions, or nonlinear relationships that the fitted model does not represent?

## Current Core Files

- `misfit.tex`: main LaTeX manuscript.
- `SEMmisfit.R`: current simulation and diagnostic script.
- `simtable.tex`: generated simulation table.
- `sim.pdf`: generated simulation figure.
- `SEMmisfit.bib`: bibliography.
- `simresults.rda`: saved simulation results.

Do not assume the script, saved results, table, and manuscript are currently consistent. Work package WP05 owns that alignment.

## Completed Roadmap

1. Improve and narrow the misfit detection approach before hardening the package interface.
2. Design a small package API around the surviving diagnostic set.
3. Reframe the manuscript around diagnostic decomposition rather than a broad anti-SEM or anti-normality argument.
4. Replace placeholder follow-up content with a concrete detect-classify-localize-model-report workflow.
5. Rebuild simulations after the diagnostic targets and defaults are settled.
6. Strengthen literature integration.
7. Add one empirical example with interpretable residual structure.
8. Run a final integration pass for terminology, claims, citations, outputs, and compilation.

## Active Revision Roadmap

The initial roadmap produced a usable prototype and manuscript draft. The next
phase corrects the residual representations and establishes estimation-aware
calibration before the work is treated as methodologically final.

1. Define separate nulls for moment fit, full Gaussian observed-data fit, and
   variable-specific conditional adequacy.
2. Use symmetric/ZCA whitening only for global diagnostics; reserve Cholesky
   innovations for an explicitly order-dependent optional comparison.
3. Localize with variable-specific conditional residuals rather than pairs of
   whitened coordinates.
4. Implement and compare conditional permutation, fixed-parameter parametric,
   and refit parametric-bootstrap calibration.
5. Run a calibration-focused simulation with computational cost determined by
   statistical precision, not convenience.
6. Revise the manuscript and finish a new integration pass only after the
   validated implementation and outputs exist.

## Working Rules For Future Sessions

- Begin by reading `CURRENT_HANDOFF.md`.
- Then read only the relevant work package in `work_packages/`.
- Keep edits inside the active work package unless the package explicitly says otherwise.
- If a dependency is discovered, record it in `CURRENT_HANDOFF.md` instead of expanding scope without reason.
- End every session by updating `CURRENT_HANDOFF.md`.
- Do not regenerate the new simulation outputs until WP09A, WP09B, and WP10A
  are complete and WP10B is active.
- Do not rewrite the full manuscript in one pass. Work section by section.

## High-Level Risks

- Overclaiming: residual diagnostics are complementary to SEM fit, not replacements.
- Calibration: residuals are computed from fitted expectations, so naive p-values may be miscalibrated.
- Interpretation: omnibus tests can detect real deviations without telling what to change.
- Scope: package support for every SEM ecosystem is too large for v1.
- Example strength: the empirical example must be interpretable, not merely statistically significant.
