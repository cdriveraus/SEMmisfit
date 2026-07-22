# WP07: Empirical Example

## Goal

Find and develop one empirical SEM example where residual diagnostics reveal interpretable structure not adequately represented by the initial SEM.

## Primary Files

- New analysis script or vignette file.
- `misfit.tex`
- `SEMmisfit.bib`
- Package prototype files from WP02, if available.

## Read First

- `CURRENT_HANDOFF.md`
- `detection_method_notes.md`
- `work_packages/WP02_software_api.md`
- Current manuscript example/simulation sections.

## Do Not Change

- Do not choose an example only because it gives a significant p-value.
- Do not use restricted data unless permissions and reproducibility are clear.
- Do not let the empirical example drive unsupported changes to the diagnostic defaults.

## Required Output

Identify candidate datasets and select one. Prefer:

- familiar open SEM datasets,
- enough sample size for nonparametric diagnostics,
- variables with interpretable possible nonlinear or subgroup structure,
- easy reproducibility in R.

Develop a concise analysis:

- fit a conventional SEM,
- report classical fit,
- run residual diagnostics,
- localize the strongest signal,
- show a simple plot,
- fit or sketch a plausible follow-up model,
- explain what the example demonstrates and what it does not prove.

## Acceptance Criteria

- Dataset is open and reproducible.
- Example is interpretable to SEM readers.
- Diagnostics use the WP01 default set.
- The follow-up model is substantively plausible, not just mechanically optimized.
- Manuscript text and figures can be added without overwhelming the paper.

## Handoff Notes

Update `CURRENT_HANDOFF.md` with:

- datasets considered,
- selected dataset and source,
- scripts created,
- diagnostic findings,
- whether the example is strong enough for the manuscript.

