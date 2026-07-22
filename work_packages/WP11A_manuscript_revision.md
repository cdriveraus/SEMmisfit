# WP11A: Manuscript Methods And Evidence Revision

## Goal

Revise the manuscript so its theory, terminology, simulation claims, software
description, and empirical interpretation match the validated WP09--WP10
implementation.

## Primary Files

- `misfit.tex`
- `SEMmisfit.bib`
- `detection_method_notes.md`
- `CURRENT_HANDOFF.md`

## Read First

- `CURRENT_HANDOFF.md`
- Completed WP09A, WP09B, WP10A, and WP10B handoff entries.
- Full `misfit.tex`, generated simulation table/figure, and revised prototype
output documentation.

## Do Not Change

- Do not add new diagnostics or alter generated numerical outputs.
- Do not rebuild package architecture.
- Do not replace the Wage example with a new empirical dataset in this
package.

## Required Output

- Rewrite residual-method sections around the three diagnostic targets and the
distinction between symmetric-whitened global checks, conditional localization,
and optional Cholesky innovations.
- Add one central methods table specifying each residual object, null,
statistic, calibration, invariance, valid interpretation, and invalid
interpretation.
- Explain conditional permutation, fixed-parameter simulation, and refit
bootstrap without treating them as interchangeable p-values.
- Replace prior simulation claims with the WP10B results and Monte Carlo
uncertainty.
- Reframe Wage as a reproducible workflow illustration, not a confirmatory or
general empirical validation; report diagnostic effect-size change across its
follow-up models.
- Add or update citations only where they are used in the revised text.

## Acceptance Criteria

- Every reported diagnostic has an exact null and calibrated interpretation.
- No text calls all rejections nonlinear SEM misspecification.
- No text uses whitened-coordinate pairs as substantively identified original
variable pairs.
- All simulation values and captions match WP10B artifacts.
- All cited keys resolve in the bibliography.

## Handoff Notes

Record changed sections, claims deliberately removed or softened, citations
added, compile status if run, and the exact WP11B checks still required.
