# WP11B: Final Integration And Reproducibility

## Goal

Perform the final consistency, reproducibility, and compilation pass after the
revised implementation, validation simulation, and manuscript are complete.

## Primary Files

- `misfit.tex`
- `SEMmisfit.bib`
- `SEMmisfit.R`
- `sem_misfit_api_prototype.R`
- `simtable.tex`
- `sim.pdf`
- `CURRENT_HANDOFF.md`

## Read First

- `CURRENT_HANDOFF.md`
- `PROJECT_PLAN.md`
- All completed WP09--WP11A handoff entries.
- Full manuscript and the latest generated outputs.

## Do Not Change

- Do not add new methods, calibration modes, simulation conditions, or
empirical examples.
- Do not silently rerun or replace outputs if their provenance is unclear.

## Required Output

- Check terminology and numerical consistency across code, manuscript,
simulation artifacts, figures, captions, and handoff.
- Run the documented prototype tests and verify that the simulation can be
regenerated from its recorded command/environment.
- Compile with `latexmk -pdf -interaction=nonstopmode -halt-on-error
misfit.tex` and resolve fatal errors, citations, references, and overfull-box
regressions caused by the revision.
- Update the reproducibility note with exact scripts, calibration modes,
replication counts, and package scope.
- Record remaining scientific limitations separately from implementation
blockers.

## Acceptance Criteria

- Manuscript compiles with no fatal errors, unresolved citations, unresolved
references, or new overfull boxes.
- All reported outputs are traceable to scripts and recorded settings.
- Handoff distinguishes submission-readiness from outstanding future package
work.
- `CURRENT_HANDOFF.md` names no unresolved blocker required for this revised
manuscript state.

## Handoff Notes

Record build command, test commands, warnings, changed files, provenance
checks, residual risks, and the recommended next task or review stage.
