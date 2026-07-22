# Work Packages

This directory splits the project into continuation-friendly units. Each future LLM session should read `../CURRENT_HANDOFF.md`, then the one work package it is executing.

## Completed Package Order

1. `WP01_detection_method.md`: improve and narrow the diagnostic approach.
2. `WP02_software_api.md`: design and prototype the R package interface.
3. `WP03_manuscript_reframe.md`: revise title, abstract, introduction, and core argument.
4. `WP04_placeholder_section.md`: replace the AI-generated follow-up section.
5. `WP05_simulation_rebuild.md`: align simulation code, outputs, and manuscript claims.
6. `WP06_literature_links.md`: strengthen literature integration.
7. `WP07_empirical_example.md`: find and develop an empirical example.
8. `WP08_final_integration.md`: final manuscript consistency and compile pass.

## Revision Package Order

The completed WP01--WP08 sequence established the current prototype and
manuscript. The following packages supersede its former final-integration
state. Run them in order and keep each chat limited to one package.

9. `WP09A_residual_targets.md`: define diagnostic nulls and replace default Cholesky whitening with symmetric whitening for global checks.
10. `WP09B_conditional_localization.md`: implement conditional-residual localization and representation-level tests.
11. `WP10A_bootstrap_calibration.md`: implement fixed-parameter and refit parametric calibration.
12. `WP10B_calibration_simulation.md`: rebuild the validation simulation and generated outputs.
13. `WP11A_manuscript_revision.md`: revise manuscript methods, results, and terminology from the validated implementation.
14. `WP11B_final_integration.md`: perform final consistency, reproducibility, and compilation checks.

## Dependency Notes

- WP01 should happen before WP02 and WP05.
- WP02 can begin as an API sketch before WP01 is final, but package defaults should not be frozen until WP01 is complete.
- WP03 and WP04 can proceed after WP01 has at least a provisional diagnostic taxonomy.
- WP06 can proceed in parallel with WP03 or WP04 if it only adds notes and bibliography candidates.
- WP07 can start early as a search task, but the final example should use the WP01 diagnostic set.
- WP08 must be last.
- WP09A precedes all new implementation work because it defines the residual
  objects and their targets.
- WP09B depends on WP09A.
- WP10A depends on WP09A and WP09B; WP10B depends on WP10A.
- WP11A depends on the validated WP10B outputs; WP11B is last.

## Standard Work Package Structure

Each file uses this structure:

- Goal
- Primary files
- Read first
- Do not change
- Required output
- Acceptance criteria
- Handoff notes

Do not broaden a work package just because adjacent work is visible. Record adjacent needs in `../CURRENT_HANDOFF.md`.
