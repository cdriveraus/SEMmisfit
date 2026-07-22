# WP07 Nonlinear Candidate Search

Date: 2026-07-13.

Scope: empirical-example search focused on nonlinear residual dependence across variables, while still recording distributional, row-distance, and marginal violations separately.

## Why The Previous Screen Was Coarse

The first WP07 screen used `B = 19` resamples for dependence tests. With that setting, the smallest attainable nonzero permutation p-value is about `.05`, so results could rank statistics but could not provide useful evidence strength or adjusted pairwise localization. It was a triage pass, not a candidate-selection pass.

## Stronger Search Setup

- Script: `empirical_example_nonlinear_search.R`.
- Full complete-case data were used; no row subsampling.
- Global residual independence: permutation dHSIC with `B = 199`.
- Pairwise localization: all pairs ranked by distance correlation on full data; the top five pairs per candidate were tested with `energy::dcov.test(..., R = 199)` and Holm-adjusted within the inspected top-pair set.
- Other violations recorded separately: energy multivariate-normal residual distribution, row-distance Anderson--Darling, and marginal Anderson--Darling checks.
- Visual follow-up script: `empirical_example_candidate_plots.R`.

Generated outputs:

- `empirical_example_nonlinear_candidates.csv`
- `empirical_example_nonlinear_top_pairs.csv`
- `empirical_example_nonlinear_marginal.csv`
- `empirical_example_nonlinear_smooth_checks.csv`
- `empirical_example_figures/wage_logwage_residual_by_age.png`
- `empirical_example_figures/nhanes_bpsys_residual_by_age.png`
- `empirical_example_figures/nhanes_bpsys_residual_by_bmi.png`

## Candidate Results

| Candidate | Source | N | Model type | Classical fit | Nonlinear residual-dependence evidence | Other violations | Judgment |
|---|---:|---:|---|---|---|---|---|
| WageLogPath | `ISLR2::Wage` | 3000 | Observed-variable path model: `logwage ~ age + year + education_num`; predictors freely covary | Saturated mean/covariance model: df = 0, CFI = 1, RMSEA = 0, SRMR near 0 | Global dHSIC p = .005. Top pairs: logwage--education_num, logwage--age, logwage--year all p = .005. GAM check: nonlinear age smooth edf = 4.42, p < machine display precision. | Strong global distribution, row-distance, and marginal signals. | Best current example. It cleanly shows that a linear SEM/path model can reproduce mean/covariance structure while residual diagnostics reveal nonlinear cross-variable dependence. |
| NHANESBpPath | `NHANES::NHANES`, adults | 6794 | Observed-variable health path model: systolic BP regressed on age, BMI, diastolic BP, cholesterol, pulse; predictors freely covary | Saturated mean/covariance model: df = 0, CFI = 1, RMSEA = 0, SRMR near 0 | Global dHSIC p = .005. Top pairs include Age--BPDiaAve, BPSysAve--Age, BPSysAve--BMI all p = .005. GAM checks: systolic BP age smooth edf = 4.56, BMI smooth edf = 1.01, both p < machine display precision. | Strong distributional and marginal signals; clinical variables have heavier tails/outliers. | Strong second candidate. It is more substantively health-relevant but the signal mixes nonlinear age/BP structure with outliers and physiological heterogeneity. |
| MSQAffectCfa | `psychTools::msq` | 1997 | Four-factor affect-item CFA | CFI = .958, RMSEA = .062, SRMR = .052 | Global dHSIC p = .005. Top pairs are near-synonym affect items, e.g., nervous--scared and active--energetic, all p = .005. | Strong distributional, row-distance, and marginal signals. | Useful for local item dependence, but not the best nonlinear example because it mostly shows redundant item content and floor/skew effects. |
| SPIBigFiveCfa | `psychTools::spi` | 4000 | Five-factor personality-item CFA | CFI = .961, RMSEA = .050, SRMR = .043 | Global dHSIC p = .005. Top pairs are same-facet item pairs, e.g., q_90--q_1763 and q_4252--q_1989, all p = .005. | Strong distributional, row-distance, and marginal signals. | Technically strong, but weaker narratively because item labels require dictionary mapping and the signal likely reflects item local dependence. |

## Recommendation

Use `ISLR2::Wage` as the leading empirical example unless the paper must stay strictly inside psychological/latent-variable data. It gives the cleanest conceptual demonstration:

1. Fit an ordinary linear observed-variable path model.
2. Classical SEM mean/covariance fit is perfect because the model is saturated for the chosen observed variables.
3. Residual-dependence diagnostics still reject, with localized dependence involving log wage and age/education.
4. A simple residual plot shows the substantive pattern: log-wage residuals rise into midlife and fall later, indicating a nonlinear age effect that the linear path model did not represent.
5. A follow-up model with a spline or quadratic age term is substantively plausible and directly matched to the diagnostic pattern.

Use `NHANES::NHANES` as the second-best option if a health/biometric example is preferable. It is larger and open, and the age/BP pattern is visually clear, but distributional/outlier and physiological-heterogeneity issues are more entangled.

Do not prioritize MSQ/SPI for the main paper example unless the desired story is local item dependence rather than nonlinear cross-variable dependence.
