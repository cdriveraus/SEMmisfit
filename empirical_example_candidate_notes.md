# WP07 Empirical Example Candidate Screen

Date: 2026-07-10.

Scope: candidate selection only. No manuscript example has been selected yet.

## Screening Setup

- Script: `empirical_example_screen.R`.
- Diagnostics: WP01 default set through `sem_misfit_api_prototype.R`.
- SEM engine: `lavaan::cfa()` with ML, mean structure, `std.lv = TRUE`, listwise complete data.
- Resampling: `B = 19` for the quick screen. This keeps runtime manageable but gives coarse dependence-test p-values. Treat p-values as triage only.
- No row subsampling was used. All complete cases for each candidate were analyzed.

Generated outputs:

- `empirical_example_candidates.csv`
- `empirical_example_candidate_tests.csv`
- `empirical_example_candidate_top_pairs.csv`
- `empirical_example_candidate_top_variables.csv`

## Candidate Results

| Candidate | Source | N analyzed | Variables | Classical fit | Residual signal | Initial judgment |
|---|---:|---:|---:|---|---|---|
| BFI | `psychTools::bfi` | 2499 | 20 | CFI = .837, RMSEA = .076, SRMR = .067 | Strong global distribution, row-distance, and marginal signals; dHSIC screen at p = .05; top pairs include N3--N4, N1--N2, A2--A3 | Familiar and reproducible, but the example may mostly show ordinal/Likert distribution and same-facet item local dependence rather than a substantively interesting nonlinear SEM problem. |
| SPI | `psychTools::spi` | 4000 | 15 | CFI = .961, RMSEA = .050, SRMR = .043 | Strong global distribution, row-distance, and marginal signals; dHSIC screen at p = .05; top pairs mostly within Big Five item clusters | Best classical fit among large candidates and largest N. Less familiar to SEM readers; item-code labels are less interpretable unless mapped through the dictionary. |
| MSQ | `psychTools::msq` | 1997 | 16 | CFI = .958, RMSEA = .062, SRMR = .052 | Strong global distribution, row-distance, and marginal signals; dHSIC screen at p = .05; top pairs are affect-adjective clusters such as active--energetic and nervous--scared | Substantively interpretable affect data, but strongest localization may be redundant item content and marginal skew/floor effects. |
| SAI | `psychTools::sai`, `time == 1` | 2948 | 16 | CFI = .824, RMSEA = .121, SRMR = .082 | Strong global distribution, row-distance, and marginal signals; dHSIC screen at p = .05; top pairs include worrying--worried, secure--confident, tense--upset | Large and interpretable state-anxiety data, but classical fit is already poor and floor effects likely dominate. Good for distributional diagnostics, weaker as a clean nonlinear-residual example. |
| Holzinger-Swineford | `psychTools::holzinger.swineford` | 301 | 11 | CFI = .959, RMSEA = .067, SRMR = .047 | Global distribution screen rejects, but row-distance and marginal diagnostics are not strong | Familiar SEM dataset, but too small for the intended nonparametric residual-dependence demonstration. |

## Ranking For The Paper

1. `psychTools::msq`: best balance of open data, interpretable variables, reasonable model fit, and substantive story around affect-item residual structure. Concern: the current strongest signals look like local dependence among near-synonyms and marginal skew/floor effects, not a clearly nonlinear structural omission.
2. `psychTools::spi`: strongest technical candidate by N and conventional fit. Concern: less familiar, harder item labels, and likely same-facet local dependence unless dictionary labels make the pattern substantive.
3. `psychTools::bfi`: most familiar large personality candidate. Concern: conventional fit is not good and residual signals are likely dominated by ordinal response distributions and item wording/facet local dependence.
4. `psychTools::sai`: large and easy to explain, but the initial model is too poor and floor effects are likely too dominant for the intended example.
5. Holzinger-Swineford: familiar fallback, but too small and not strong enough for nonlinear-residual diagnostics.

## Recommendation Before Developing The Example

Do not yet commit this empirical example to the manuscript. The best current options are MSQ or SPI, but neither yet demonstrates an especially interesting nonlinear SEM problem. The next verification step should inspect residual plots for the top MSQ and SPI pairs and, if needed, test a more substantively targeted model where the diagnostic pattern plausibly reflects local dependence, subgroup heterogeneity, or nonlinear association rather than only nonnormal ordinal responses.
