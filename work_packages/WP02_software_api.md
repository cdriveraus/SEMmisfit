# WP02: Software API And Package Prototype

## Goal

Design a small R package interface for residual-based SEM misfit diagnostics. The package should make the diagnostics from WP01 usable from common SEM packages without committing to broad package support too early.

## Primary Files

- New package/prototype files, if created.
- `SEMmisfit.R` as the source of current prototype logic.
- `detection_method_notes.md` from WP01.

## Read First

- `CURRENT_HANDOFF.md`
- `work_packages/WP01_detection_method.md`
- `detection_method_notes.md`, if it exists.
- In `SEMmisfit.R`, read the fitted expectation extraction and residual standardization code.

## Do Not Change

- Do not change manuscript prose.
- Do not regenerate simulation outputs.
- Do not implement readers for every SEM package. Keep v1 narrow.

## Required Output

Produce an API sketch or prototype with these concepts:

- `as_sem_expectation(object, data = NULL, ...)`
  - returns observed data, implied mean vector, implied covariance matrix, variable names, sample size, group metadata if available.
- `sem_misfit(object, data = NULL, tests = "default", calibration = "auto", ...)`
  - returns residuals, diagnostics, statistics, p-values or calibrated evidence, variable-pair results, and plot-ready data.
- Reader methods for v1:
  - lavaan,
  - OpenMx,
  - generic list input for tests independent of SEM software.
- Output object class:
  - print method with compact summary,
  - data-frame extractor for test results,
  - plot methods for pairwise heatmap and rowwise misfit.

## Acceptance Criteria

- API accepts fitted SEM objects and a generic expectation object.
- v1 scope is explicit: continuous observed variables, single-group first, complete-data/listwise handling.
- Defaults come from WP01.
- Readers are thin wrappers around one shared diagnostic engine.
- The prototype can reproduce the current simulation diagnostics in principle, even if full implementation is deferred.

## Handoff Notes

Update `CURRENT_HANDOFF.md` with:

- API decisions,
- files created,
- supported SEM sources,
- unsupported cases,
- next package implementation step.

