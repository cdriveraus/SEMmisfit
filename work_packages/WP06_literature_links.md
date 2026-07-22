# WP06: Literature Links

## Goal

Strengthen the manuscript's connection to existing SEM fit, residual diagnostics, nonlinear SEM, and nonparametric independence-testing literatures.

## Primary Files

- `misfit.tex`
- `SEMmisfit.bib`
- Optional new notes file: `literature_notes.md`

## Read First

- `CURRENT_HANDOFF.md`
- `misfit.tex`
- `SEMmisfit.bib`
- `detection_method_notes.md`, if it exists.

## Do Not Change

- Do not add citations that are not used in the manuscript or notes.
- Do not do a broad citation dump.
- Do not rewrite entire manuscript sections unless the active task explicitly requires it.

## Required Output

Create or update `literature_notes.md` with:

- Existing bibliography entries that should be used more strongly.
- Missing references to add, grouped by purpose:
  - classical SEM fit and residual diagnostics,
  - robust/non-normal SEM,
  - HFI and nonlinear SEM,
  - dHSIC/HSIC, distance covariance, energy tests, HHG if useful,
  - SEM trees, LOSEM, model-based recursive partitioning,
  - Bayesian posterior predictive checks,
  - empirical example source if identified.
- Specific manuscript locations where citations should be inserted.

Patch `SEMmisfit.bib` only for references that are clearly needed.

## Acceptance Criteria

- Literature notes identify which claims each citation supports.
- Manuscript citation gaps are concrete, not generic.
- Added bibliography entries compile under biblatex/APA style.
- The paper's novelty is positioned as integration and software workflow, not invention of all diagnostics.

## Handoff Notes

Update `CURRENT_HANDOFF.md` with:

- citations added,
- citation gaps remaining,
- papers that should be read before finalizing claims,
- manuscript sections most affected.

