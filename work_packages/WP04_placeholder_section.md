# WP04: Replace Placeholder Follow-Up Section

## Goal

Replace the AI-generated content in `\section{Non-linear Misspecification Found -- What Now?}` with a concrete workflow for what researchers should do after residual diagnostics indicate misfit.

## Primary Files

- `misfit.tex`
- `SEMmisfit.bib`
- `detection_method_notes.md`

## Read First

- `CURRENT_HANDOFF.md`
- `detection_method_notes.md`, if it exists.
- In `misfit.tex`, read `Detecting Nonlinear Misfit in SEM`, `Simulation Results`, `Non-linear Misspecification Found -- What Now?`, and `Conclusion`.

## Do Not Change

- Do not rewrite the introduction unless needed for a cross-reference.
- Do not alter simulation code or outputs.
- Do not present model modification as automatic repair.

## Required Output

Replace the placeholder section with a workflow:

1. Detect: run diagnostics and identify which null failed.
2. Classify: separate marginal distribution, joint distribution, residual dependence, heterogeneity, and covariance/mean misfit.
3. Localize: use pairwise heatmaps, variable-wise checks, rowwise scores, moderator scans, and residual plots.
4. Model follow-up: consider theory-guided nonlinear paths, interactions, splines, moderators, mixtures/classes, SEM trees, LOSEM, Bayesian checks, or robust distributional models.
5. Report: describe diagnostics as sensitivity checks and exploratory guides, not automatic proof of the correct alternative model.

Tie follow-up options to existing or newly added literature where possible.

## Acceptance Criteria

- The section is no longer generic.
- Every recommended follow-up is connected to the diagnostic signal it addresses.
- The section cautions against data-dredging and automatic model repair.
- The section gives a reader a practical path from detected misfit to interpretable next steps.

## Handoff Notes

Update `CURRENT_HANDOFF.md` with:

- section replaced,
- citations added or still needed,
- unresolved conceptual concerns,
- whether the conclusion now needs revision.

