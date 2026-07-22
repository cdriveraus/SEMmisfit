# Current Handoff

## Status

Date initialized: 2026-07-09.
Latest update: 2026-07-21, WP09B completed. WP01--WP08 remain a completed historical baseline, but their prior final-integration state is superseded by the new residual-representation and calibration work. Begin with WP10A only.

This repository contains a manuscript on nonlinear SEM misfit diagnostics plus prototype analysis/package artifacts. The coordination scaffold has been created so future LLM sessions can work in isolated work packages.

WP01 added a methodological notes artifact. WP02 added a prototype API file. WP03 reframed manuscript prose around residual diagnostic decomposition. WP04 replaced the follow-up/model-revision placeholder section. WP05 rebuilt the simulation workflow and regenerated `simresults.rda`, `simtable.tex`, and `sim.pdf` as illustrative outputs aligned with the WP01 diagnostic taxonomy. WP06 added a focused literature map and patched targeted citations into the manuscript. WP07 selected and developed the Wage empirical example. WP08 completed the final integration and compile pass.

## 2026-07-20 Wage Software Workflow Revision

- Retained the Wage dataset and did not add an SEM-tree dependency. The manuscript explains that the smooth age pattern is not evidence for discrete moderator-defined subgroups.
- Added a software implementation subsection with executable `lavaan` -> `sem_misfit()` -> accessor/plot code and explicit prototype scope/calibration limits.
- Added `sem_misfit_variables()` and `plot(..., type = "conditional_residual", variable = ..., against = ...)` without changing the expectation interface or diagnostic definitions.
- Replaced the spline-led main figure with a package-generated pairwise heatmap and conditional-residual plot.
- Added matched conditional `lavaan` models. Freeing the quadratic age path gave $\Delta\chi^2(1)=196.94$, `p < .001`, and changed AIC from 2091.19 to 1896.25.
- Age--residual distance correlation fell from .129 to .063 in the quadratic SEM but remained detectable (`p = .001`); the spline sensitivity result was .048 (`p = .022`). The manuscript therefore calls the modification substantial but incomplete.
- Added a two-outcome `lavaan` syntax example to show how observed basis terms generalize to multivariate path systems, while deferring case-specific conditional expectations as future API work.
- New artifacts: `wage_empirical_cases.csv`, `wage_empirical_sem_comparison.csv`, `wage_empirical_diagnostics.rds`, `empirical_example_figures/wage_pairwise_heatmap.png`, `empirical_example_figures/wage_conditional_residual_by_age.png`, and `empirical_example_figures/wage_quadratic_residual_by_age.png`.
- Validation: full `Rscript empirical_example_wage.R` completed; targeted non-diagonal whitening/accessor/plot/follow-up checks passed; final PDF pages 9--10 were visually inspected after rendering.

## Latest Decisions

- Treat the current fixed-fit residual permutation approach as conditional
  exploratory calibration, not unconditional SEM-level inference.
- Retain that conditional approach for comparison while adding
  fixed-parameter and refit parametric-bootstrap calibration.
- Use refit bootstrap as the preferred inferential standard when implemented;
  simulation cost is acceptable when needed to establish calibration.
- Do not use Cholesky-whitened coordinate pairs as original-variable
  localization. Use symmetric whitening for global checks and conditional
  residuals for variable-labelled localization.

- Frame the project as residual-based diagnostic decomposition for SEM, not as a replacement for classical SEM fit.
- Keep package v1 narrow: continuous observed variables, single-group first, complete-data/listwise diagnostics, model-implied observed means and covariance as the main expectation interface.
- Treat simulation output as provisional until the detection approach is rebuilt around the WP01 diagnostic taxonomy.
- Distinguish distributional diagnostics from residual-dependence diagnostics throughout the package and manuscript.
- Use whitened observed residuals as the primary v1 residual representation.
- Use conditional residuals mainly for variable-labeled localization and plotting, not as the main global test basis.
- Use reader methods only to create a shared `sem_expectation` object; all diagnostics should run through one shared engine.
- Keep v1 readers to generic list input, lavaan, and OpenMx.
- Keep bootstrap calibration in the API vocabulary but do not treat it as implemented until a refit-capable bootstrap engine exists.
- WP02 prototype whitening uses base R's upper Cholesky factor in row orientation: `raw %*% solve(chol(sigma))`, implemented via `backsolve(..., transpose = TRUE)`.

## Completed Work Package

Completed: `work_packages/WP01_detection_method.md`.

Output: `detection_method_notes.md`.

Chosen default diagnostics:

- Baseline SEM moment fit: standard SEM chi-square/global fit output as comparator for mean/covariance misfit.
- Global residual distribution test: energy-distance multivariate goodness-of-fit test or equivalent omnibus test of whitened residuals against standard multivariate normality; alternatively rowwise Mahalanobis-distance Anderson-Darling as a simpler scalar check.
- Global residual independence test: joint dHSIC-style test on whitened residual columns as the primary nonlinear-dependence diagnostic.
- Pairwise residual independence localization: pairwise distance covariance or pairwise HSIC/dHSIC with multiplicity adjustment.
- Marginal residual normality localization: univariate Anderson-Darling checks on whitened residual coordinates, clearly labeled as distributional localization.
- Case-level residual distance: rowwise squared Mahalanobis distances with reference quantiles; inferential p-values only if calibrated.

Diagnostics excluded from v1 defaults unless later evidence supports them:

- current HFI scalar as a p-value,
- mutual-information tests based on discretization,
- aggregate residual-rest mutual-information variants,
- Levene-style tests using residual-derived groups,
- overlap-coefficient tests against simulated references,
- naive chi-square p-values for rowwise distances as confirmatory inference.

Completed: `work_packages/WP02_software_api.md`.

Output: `sem_misfit_api_prototype.R`.

API decisions:

- `as_sem_expectation(object, data = NULL, ...)` is an S3 generic.
- `as_sem_expectation.list()` accepts generic expectation input with observed data, implied mean, implied covariance, optional variable names, source label, and optional SEM fit statistics.
- `as_sem_expectation.lavaan()` is a thin reader for single-group lavaan objects, using `lavaan::lavInspect()` for data and implied moments.
- `as_sem_expectation.MxModel()` is a thin OpenMx reader that searches fitted fit-function algebra attributes for `expMean` and `expCov`, matching the current `SEMmisfit.R` prototype.
- `sem_misfit(object, data = NULL, tests = "default", calibration = "auto", ...)` coerces to a `sem_expectation`, computes raw, whitened, conditional, and row-distance residuals, then returns one `sem_misfit` object.
- Output includes global test rows, variable-level rows, pairwise rows, case-level row distances, residual matrices, and plot-ready data.
- S3 helpers include `print.sem_misfit()`, `as.data.frame.sem_misfit()`, `plot.sem_misfit()`, `sem_misfit_pairs()`, and `sem_misfit_cases()`.
- Critical correction after WP02: `compute_sem_residuals()` now whitens residuals in the correct Cholesky direction for non-diagonal covariance, so row distances match `mahalanobis(data, mean, covariance)`.

Supported SEM sources for v1:

- generic list expectation objects,
- lavaan single-group continuous observed-variable fits,
- OpenMx single-group continuous observed-variable fits with accessible fitted expectation attributes.

Unsupported cases for v1:

- multiple groups,
- ordinal/categorical indicators,
- missing-data diagnostics beyond listwise/complete data,
- mixture models and multilevel SEM,
- robust corrections and estimator-specific residual nulls,
- automated readers for every SEM package,
- fully refit parametric bootstrap calibration.

Completed: `work_packages/WP03_manuscript_reframe.md`.

Output: patched `misfit.tex`.

Sections changed:

- title and short title,
- abstract,
- introduction,
- `SEM & Model Fit`,
- `Residuals in SEM`,
- `Detecting Nonlinear Misfit in SEM`,
- simulation-results interpretation paragraph,
- conclusion interpretation paragraph.

Claims softened:

- Replaced the broad "Great SEM Delusion" / inappropriate-reference framing with residual diagnostics for SEM misfit beyond mean and covariance fit.
- Classical SEM fit is now described as useful for its intended target: mean and covariance reproduction.
- Normality violations are no longer treated as direct evidence of nonlinear SEM misspecification.
- Multivariate normality tests and HFI-style summaries are framed as broad distributional or residual-variance diagnostics, not source-specific nonlinear-dependence tests.
- Residual-dependence diagnostics such as distance covariance/dHSIC are separated from residual distribution checks.
- Localization diagnostics are framed as triage for follow-up modeling, not automatic model revision advice.

Citations needed or to review later:

- The revised text uses existing citations for SEM fit, nonnormal psychological data, dHSIC, energy tests, and independence-test review.
- Future literature pass should check whether additional citations are needed for residual whitening/Mahalanobis diagnostics, parametric-bootstrap calibration of fitted residual diagnostics, distance covariance specifically, and tree/moderated/flexible SEM follow-up approaches.

Next manuscript section to revise:

- `Simulation Study - Misfit Detection` should be rebuilt after WP05 regenerates the simulation design/results around the WP01 diagnostic taxonomy. Until then, current simulation claims should remain provisional.

Completed: `work_packages/WP04_placeholder_section.md`.

Output: patched `misfit.tex`.

Section replaced:

- `Non-linear Misspecification Found -- What Now?`

WP04 changes:

- Replaced the generic AI-generated follow-up section with a diagnostic workflow: detect the failed residual null, classify the diagnostic signal, localize the pattern, choose follow-up models matched to the signal, and report diagnostics as sensitivity evidence.
- Connected follow-up options to diagnostic signals: mean/covariance misfit to ordinary SEM revisions; marginal distributional misfit to transformations, heavy-tailed/robust distributions, outlier modeling, and sensitivity analyses; residual dependence to nonlinear paths, latent interactions, product terms, or splines; moderator-varying residual patterns to moderated SEM and LOSEM; subgroup/case-distance patterns to mixtures, latent classes, and SEM trees; complex approximate restrictions to Bayesian SEM/posterior predictive checking.
- Added citations already present in `SEMmisfit.bib`: `rizzo2016energy`, `pfister2018kernelbased`, `karch2024pearsons`, `briley2015nonparametric`, `brandmaier2013structural`, and `muthen2012bayesian`.
- Revised the conclusion's follow-up paragraph to match the new section's cautionary framing.

Citations still needed or to review later after WP04:

- Robust/heavy-tailed distributional SEM or robust residual modeling.
- Latent interaction/product-term SEM and spline/nonlinear SEM implementations.
- Mixture/latent-class SEM as follow-up for subgroup or row-distance patterns.
- Parametric-bootstrap or posterior-predictive calibration references specific to fitted SEM residual diagnostics.

Unresolved conceptual concerns after WP04:

- The follow-up section now gives a practical diagnostic-to-modeling path, but the exact default calibration language should be revisited after WP05 tests naive, permutation, and parametric-bootstrap behavior.
- The distinction between whitened residual-coordinate localization and original-variable interpretation may need a short example after the software API stabilizes.
- Conclusion no longer needs immediate WP04-specific revision, but it should be revisited after WP05 rebuilds the simulation section and results.

Completed: `work_packages/WP05_simulation_rebuild.md`.

Output: patched `SEMmisfit.R`, patched simulation section in `misfit.tex`, regenerated `simresults.rda`, `simtable.tex`, and `sim.pdf`.

Command used:

- `Rscript SEMmisfit.R`

Runtime and package issues:

- First run stopped because R had no CRAN mirror configured; `SEMmisfit.R` now sets `options(repos = c(CRAN = "https://cloud.r-project.org"))`.
- R emitted package-build warnings for `energy` and `dHSIC` being built under R 4.6.1, but the run completed.
- A PowerShell stdin pipeline to `Rscript -` hit a BOM/encoding parse error; direct `Rscript -e` and `Rscript SEMmisfit.R` worked.

WP05 changes:

- Replaced the old OpenMx/parallel/HFI/MI simulation script with a compact observed-moment simulation workflow.
- Simulation conditions now map to WP01 targets: no misfit, omitted covariance, omitted interaction, nonlinear curve, skewed marginal, heavy tails, group mean, and group variance.
- Diagnostics now match the WP01 default taxonomy: SEM chi-square moment-fit comparator, row-distance Anderson--Darling, marginal Anderson--Darling, energy multivariate-normal residual distribution test, global dHSIC residual-independence test, and pairwise distance-covariance localization.
- HFI-style values and mutual-information diagnostics were removed from generated outputs.
- Detection rates are summarized as `p <= .05`.
- `misfit.tex` now describes the simulation as an illustrative diagnostic decomposition rather than a definitive benchmark.

Key result patterns from the regenerated outputs:

- SEM chi-square detected omitted covariance most clearly, with strongest detection at `n = 500`, and showed weaker, less systematic sensitivity to nonlinear, distributional, and subgroup conditions.
- In no-misfit cells, row-distance AD and marginal AD were conservative; SEM rejected in 5.0--8.6% of replications; energy MVN rejected in 3.6--5.6%; global dHSIC rejected in 2.4--5.4%; pairwise dCov rejected in 0--2.2%.
- Energy MVN was broadly sensitive to omitted interaction, nonlinear curve, skewed marginal, heavy tails, and group variance, especially at `n = 500`.
- Marginal AD mainly responded to coordinate-focused distributional departures at larger sample size.
- Row-distance AD mainly responded to case-distance distribution changes, especially interaction, nonlinear curve, heavy tails, and group variance at larger sample size.
- Global dHSIC was sensitive to omitted interactions and nonlinear curves, especially with four variables and larger sample size; it also reacted to the four-variable omitted-covariance condition because the incorrect covariance left residual dependence after whitening.
- Pairwise dCov strongly localized the four-variable omitted-covariance and nonlinear-curve conditions, but showed little sensitivity in most eight-variable cells after multiplicity adjustment.

Remaining WP05 limitations:

- The regenerated simulation now uses 500 replications per condition and is suitable as current discussion evidence, but it is still not a final package validation study.
- Global dHSIC must not use gamma calibration as the simulation default. A 500-replication no-misfit diagnostic showed gamma calibration was badly anti-conservative in the 8-variable, n = 100 null cell; permutation calibration brought global dHSIC false detection back near nominal.
- Pairwise dCov uses permutation p-values with Bonferroni adjustment; the eight-variable conditions show that multiplicity and whitening-coordinate interpretation need more work before making localization claims.
- A future simulation should use more replications, stronger permutation settings, and ideally refit parametric bootstrap calibration before making comparative claims.

WP05 calibration rerun on 2026-07-10:

- Patched `SEMmisfit.R` to use `niter = 500`, `nperm = 199`, Windows `parallel::parLapplyLB()`, and permutation-calibrated global dHSIC.
- Regenerated `simresults.rda`, `simtable.tex`, and `sim.pdf` with 16,000 simulated datasets.
- Added `simsummary_n500.csv` as a compact full-condition detection-rate table.
- Added `null_calibration_check.R`, `null_calibration_check_results.rda`, and `null_calibration_summary_n500.csv` for the focused no-misfit calibration diagnostic.
- No-misfit false detection rates from the full rerun were:
  - `n = 100`, 4 variables: SEM .064, row-distance AD .004, marginal AD .000, energy MVN .056, global dHSIC .046, pairwise dCov .016.
  - `n = 500`, 4 variables: SEM .068, row-distance AD .000, marginal AD .000, energy MVN .036, global dHSIC .024, pairwise dCov .022.
  - `n = 100`, 8 variables: SEM .086, row-distance AD .000, marginal AD .000, energy MVN .054, global dHSIC .044, pairwise dCov .000.
  - `n = 500`, 8 variables: SEM .050, row-distance AD .002, marginal AD .000, energy MVN .052, global dHSIC .054, pairwise dCov .000.
- Interpretation: the main anti-conservative residual diagnostic was global dHSIC with gamma calibration. Permutation-calibrated global dHSIC and energy MVN are now close to nominal under no misfit. Row-distance AD, marginal AD, and Bonferroni pairwise dCov are conservative when applied to fitted residuals without refit-bootstrap calibration. The observed-moment SEM comparator is mildly liberal in some smaller null cells, especially `n = 100`, 8 variables.

Completed: `work_packages/WP06_literature_links.md`.

Output: created `literature_notes.md`, patched `misfit.tex`, patched `SEMmisfit.bib`.

Bibliography entries added:

- `bollen1992bootstrapping`: Bollen-Stine bootstrap calibration for SEM fit measures.
- `satorra1994corrections`: robust corrections to SEM test statistics and standard errors under nonnormality.
- `klein2000maximum`: LMS latent interaction method for nonlinear/interaction follow-up modeling.
- `lubke2005investigating`: factor mixture models for population heterogeneity follow-up.
- `irmer2024estimating`: recent nonlinear SEM / moderation power and `powerNLSEM` software context.
- `cox2025comparing`: recent Bayesian and structural-after-measurement comparison for latent-interaction SEM with complex data.
- `jia2026obtaining`: recent robust practical fit-index work for multiply imputed nonnormal SEM data.
- `garnierVillarreal2020adapting`: Bayesian SEM fit-index adaptation.
- `arnold2021score`: recent score-guided SEM trees.
- `silvaDiaz2025evaluation`: recent SEM forests evaluation for omitted influential covariates.

Existing bibliography entries now used more directly:

- `gerhard2017fit` for HFI as a related residual-variance/nonlinear-SEM proposal.
- `mcdonald1962general` for early nonlinear factor-analysis precedent.
- `rizzo2016energy`, `pfister2018kernelbased`, and `karch2024pearsons` for energy, dHSIC, distance covariance, and HHG-style independence-test context.
- `brandmaier2013structural`, `briley2015nonparametric`, and `muthen2012bayesian` for SEM trees, LOSEM/moderated SEM, and Bayesian SEM/posterior predictive follow-up.

WP06 manuscript sections most affected:

- Introduction: novelty positioned as integration and software workflow rather than invention of all component diagnostics.
- Introduction: recent SEM-specific layer added for nonlinear SEM, latent interactions, robust fit indices, Bayesian fit indices, SEM trees, and SEM forests.
- `SEM & Model Fit`: robust corrections and bootstrap calibration cited while preserving the mean/covariance target of classical fit indices.
- `Residuals in SEM`: rowwise residual-distance inference tied to bootstrap calibration caveat.
- `Detecting Nonlinear Misfit in SEM`: distance/kernel tests and HFI paragraphs linked to specific literature.
- `Simulation Study - Misfit Detection`: stronger calibration caveat cited.
- `Non-linear Misspecification Found -- What Now?`: follow-up model families tied to robust SEM, recent robust fit work, latent interaction, mixture, SEM tree/forest, LOSEM, and Bayesian SEM references.

Post-WP06 formatting pass:

- `misfit.tex` math delimiters were normalized from `\(...\)` and `\[...\]` to `$...$` and `$$...$$` at the user's request.

WP06 citation gaps remaining:

- Read a focused heavy-tailed/t-distribution SEM source before making detailed claims about distributional SEM follow-up; current text only makes broad robust-correction and robust-fit-index claims.
- Add a direct Mahalanobis/outlier-diagnostics reference only if rowwise distances become a major method rather than a short residual summary.
- Add original distance covariance, HSIC, or HHG references only if the manuscript gives technical method definitions rather than workflow-level citations.
- Add empirical-example citations after WP07 selects the example data source.

Papers or areas to read before finalizing claims:

- Yuan/Bentler-style robust mean and covariance structure analysis and recent robust SEM test-statistic comparisons.
- Heavy-tailed or t-distribution SEM if robust/heavy-tail follow-up becomes a concrete recommendation.
- General model-based recursive partitioning if the follow-up section expands beyond SEM trees.
- Broader posterior predictive checking references if Bayesian calibration becomes central.

## WP07 Candidate Screen

Status: candidate screening completed; empirical example not yet selected.

Files created:

- `empirical_example_screen.R`
- `empirical_example_candidate_notes.md`
- `empirical_example_candidates.csv`
- `empirical_example_candidate_tests.csv`
- `empirical_example_candidate_top_pairs.csv`
- `empirical_example_candidate_top_variables.csv`

Screening setup:

- Used `sem_misfit_api_prototype.R` with the WP01 default diagnostics.
- Fit simple `lavaan::cfa()` models with ML, mean structure, `std.lv = TRUE`, and complete/listwise data.
- Used all complete cases for each candidate; no row subsampling.
- Used `B = 19` for quick dependence-test triage, so dependence-test p-values are coarse and should not be treated as final calibration.

Datasets considered:

- `psychTools::bfi`: 2499 complete cases, 20 selected Big Five items. Familiar and open; strong residual distribution and marginal signals, but likely dominated by Likert response shape and same-facet item local dependence.
- `psychTools::spi`: 4000 complete cases, 15 selected SAPA personality items. Best technical candidate by N and conventional fit (CFI about .961, RMSEA about .050, SRMR about .043), but less familiar and item labels need dictionary mapping.
- `psychTools::msq`: 1997 complete cases, 16 affect items. Strong substantive interpretability and reasonable fit (CFI about .958, RMSEA about .062, SRMR about .052), but top localized pairs are near-synonym affect items and may mainly show local dependence plus marginal skew/floor effects.
- `psychTools::sai`, restricted to `time == 1`: 2948 complete cases, 16 state-anxiety items. Large and interpretable, but the initial two-factor model fit is poor (CFI about .824, RMSEA about .121), and floor effects likely dominate.
- `psychTools::holzinger.swineford`: 301 complete cases, 11 cognitive indicators. Familiar SEM example, but too small and not strong enough for the intended nonparametric residual-dependence demonstration.

Initial ranking:

1. `psychTools::msq`
2. `psychTools::spi`
3. `psychTools::bfi`
4. `psychTools::sai`
5. `psychTools::holzinger.swineford`

Current WP07 judgment:

- Do not yet commit an empirical example to the manuscript.
- MSQ and SPI are the best current candidates, but neither yet clearly demonstrates an especially interesting nonlinear SEM problem.
- The next verification step should inspect residual plots for the top MSQ and SPI pairs and decide whether the signal is substantively interesting enough, or whether a more targeted open dataset/problem is needed.

## WP07 Nonlinear Candidate Search

Status: stronger nonlinear-dependence search completed; empirical example still awaits user verification before manuscript integration.

Files created:

- `empirical_example_nonlinear_search.R`
- `empirical_example_nonlinear_notes.md`
- `empirical_example_candidate_plots.R`
- `empirical_example_nonlinear_candidates.csv`
- `empirical_example_nonlinear_top_pairs.csv`
- `empirical_example_nonlinear_marginal.csv`
- `empirical_example_nonlinear_smooth_checks.csv`
- `empirical_example_figures/wage_logwage_residual_by_age.png`
- `empirical_example_figures/nhanes_bpsys_residual_by_age.png`
- `empirical_example_figures/nhanes_bpsys_residual_by_bmi.png`

Search setup:

- Used full complete-case data with no row subsampling.
- Used permutation-calibrated global dHSIC with `B = 199`.
- Ranked all pairwise residual-coordinate dependencies by distance correlation, then tested the top five pairs per candidate with `energy::dcov.test(..., R = 199)` and Holm adjustment within the inspected top-pair set.
- Recorded residual-distribution, row-distance, and marginal-normality violations separately from nonlinear residual-dependence diagnostics.
- Added simple GAM smooth checks and residual plots for the two leading observed-variable path-model candidates.

Datasets/models considered:

- `ISLR2::Wage`: 3000 rows. Observed-variable path model `logwage ~ age + year + education_num`, with predictors freely covarying. This model is saturated for the selected observed variables, so classical SEM mean/covariance fit is perfect (`df = 0`, CFI = 1, RMSEA = 0, SRMR near 0). Residual diagnostics still show global residual independence failure (`p = .005`) and top pairwise dependencies involving `logwage` with `education_num`, `age`, and `year` (`p = .005` for each top tested pair). A GAM check of `logwage ~ s(age) + year + education_num` found a nonlinear age smooth (`edf = 4.42`, displayed `p = 0`). The residual plot shows log-wage residuals rising into midlife and falling later.
- `NHANES::NHANES` adults: 6794 complete cases. Observed-variable path model with systolic blood pressure regressed on age, BMI, diastolic blood pressure, direct cholesterol, total cholesterol, and pulse, with predictors freely covarying. Also saturated for mean/covariance (`df = 0`, CFI = 1, RMSEA = 0, SRMR near 0). Residual diagnostics show global residual independence failure (`p = .005`) and top dependencies involving age, blood pressure, BMI, and cholesterol. GAM checks show nonlinear age/BP structure (`edf = 4.56`, displayed `p = 0`) and a weaker BMI pattern.
- `psychTools::msq`: 1997 rows. Four-factor affect CFA still shows global dHSIC failure, but top pairs are near-synonym affect items such as `nervous--scared` and `active--energetic`; better interpreted as local item dependence plus marginal/floor effects.
- `psychTools::spi`: 4000 rows. Five-factor personality CFA still shows global dHSIC failure, but top pairs are same-facet item pairs with opaque item-code labels; better as local item dependence than nonlinear cross-variable demonstration.

Current recommendation:

- Use `ISLR2::Wage` as the leading empirical example unless the paper must stay strictly inside psychological latent-variable data. It is the cleanest demonstration that a linear SEM/path model can reproduce mean/covariance structure while residual diagnostics reveal nonlinear cross-variable dependence.
- Use `NHANES::NHANES` as the second-best option if a health/biometric example is preferable. It is larger and visually clear, but nonlinear dependence is more entangled with outliers, heterogeneity, and distributional violations.
- Do not prioritize MSQ/SPI for the main example unless the desired story is local item dependence rather than nonlinear cross-variable dependence.

## WP07 HBN Feasibility Check

Status: checked HBN as a possible p-factor psychopathology example.

File created:

- `empirical_example_hbn_feasibility.md`

Findings:

- HBN is scientifically attractive for the paper because it is transdiagnostic pediatric mental-health data and includes CBCL, which is central to youth p-factor/internalizing/externalizing work.
- Public HBN-EEG/OpenNeuro releases expose derived `p_factor`, `attention`, `internalizing`, and `externalizing` scores in `participants.tsv` across 3155 checked participants, but not raw CBCL item or syndrome-scale tables.
- The HBN Release 11 public data dictionaries confirm that raw CBCL variables exist (`CBCL_01`, `CBCL_02`, etc. for age 6-18; `CBCL_Pre_01`, etc. for preschool), but full phenotypic data require HBN DUA approval and LORIS access.
- Derived p/internalizing/externalizing scores alone are not sufficient for the manuscript question because they presuppose the factor model. To test whether apparent internalizing/externalizing factors arise from nonlinear residual dependence, we need raw CBCL indicators.

Current HBN recommendation:

- If DUA-gated data are acceptable, HBN is the best substantive next candidate: request/use LORIS access, download CBCL, and compare one-factor, two-factor, bifactor, and nonlinear-residual follow-up models.
- If the empirical example must be fully reproducible without access steps, HBN should not replace the fully open Wage example.

## WP07 Wage Empirical Example

Status: selected and developed as the current empirical example.

Files created or updated:

- `empirical_example_wage.R`
- `wage_empirical_summary.csv`
- `wage_empirical_diagnostics.csv`
- `wage_empirical_pairwise.csv`
- `wage_empirical_marginal.csv`
- `wage_empirical_cases.csv`
- `wage_empirical_followup.csv`
- `wage_empirical_sem_comparison.csv`
- `wage_empirical_smooth_checks.csv`
- `wage_empirical_diagnostics.rds`
- `empirical_example_figures/wage_linear_residual_by_age.png`
- `empirical_example_figures/wage_pairwise_heatmap.png`
- `empirical_example_figures/wage_conditional_residual_by_age.png`
- `empirical_example_figures/wage_quadratic_residual_by_age.png`
- `empirical_example_figures/wage_spline_residual_by_age.png`
- `misfit.tex`
- `SEMmisfit.bib`
- `misfit.pdf`

Selected dataset:

- `ISLR2::Wage`, 3000 observations, cited via `james2021introduction`.

Baseline model:

- Observed-variable path model on `age_decade`, centered survey year (`year_c`), and `logwage`.
- Model syntax: `logwage ~ age_decade + year_c; age_decade ~~ year_c`.
- The model is saturated for the selected mean/covariance structure (`df = 0`, CFI = 1, RMSEA = 0, SRMR near 0), which makes it a clean demonstration that classical covariance fit can be silent while residual structure remains.

Diagnostic findings:

- Global residual independence rejected with permutation dHSIC (`p = .001`, `B = 999`).
- Pairwise distance-covariance localization identified `logwage--age_decade` as the strongest pair (Holm-adjusted `p = .003`), with `logwage--year_c` also adjusted `p = .003`.
- Distributional diagnostics also rejected: global residual distribution, row-distance AD, and marginal AD all flagged the fitted Gaussian residual reference. Manuscript text states this is not a pure nonlinear-dependence example and separates these signals.

Follow-up findings:

- The linear model residuals show clear curvature over age: residual log wages rise into midlife and fall later.
- A spline follow-up `logwage ~ s(age, k = 6) + year_c` confirmed nonlinear age structure (`edf = 4.60`, `p < .001`) and reduced age-residual distance correlation from .129 to .048.
- Including education as a covariate in the spline follow-up reduced the age-residual distance-covariance test to non-significance (`p = .128`), supporting the manuscript framing that follow-up models are sensitivity evidence, not automatic repairs.

Build check:

- `latexmk -pdf -interaction=nonstopmode -halt-on-error misfit.tex` completed successfully.
- No unresolved citations or references remained in the final log scan.
- At the time of WP07, an overfull hbox warning remained at the old simulation equation line; WP08 later fixed this formatting issue.

Completed: `work_packages/WP08_final_integration.md`.

Output: patched `misfit.tex`, rebuilt `misfit.pdf`, updated `CURRENT_HANDOFF.md`.

WP08 integration changes:

- Standardized remaining manuscript wording around `nonlinear`, `modeling`, and `chi-square`.
- Repaired mojibake in the author note by using LaTeX-safe `Binzm\"uhlestrasse`.
- Added a concise reproducibility note to the author note covering simulation, empirical example, figures, and prototype diagnostic workflow scripts.
- Added a final limitations paragraph covering illustrative simulation scope, residual p-value calibration, the observed-variable Wage example, continuous single-group complete-data scope, and future package work.
- Split the group-variance simulation equation to remove the previous overfull hbox warning.

WP08 build check:

- Compile command used: `latexmk -pdf -interaction=nonstopmode -halt-on-error misfit.tex`.
- Build completed successfully and wrote `misfit.pdf`.
- Final log scan found no unresolved citations, unresolved references, fatal errors, or overfull boxes.
- Remaining warnings are ordinary underfull vbox layout warnings from LaTeX page balancing.

WP08 readiness state:

- The manuscript is internally consistent with the implemented diagnostics and generated results.
- The manuscript is suitable for coauthor/submission-readiness review, but not final package release.
- Remaining work is polish/scientific review and future package hardening, not a manuscript compilation blocker.

## Completed: `work_packages/WP09A_residual_targets.md`

Date: 2026-07-21.

Files changed:

- `sem_misfit_api_prototype.R`
- `detection_method_notes.md`
- `misfit.tex`
- `misfit.pdf`
- `CURRENT_HANDOFF.md`

WP09A changes:

- Defined exactly three diagnostic targets in the prototype metadata, notes, and manuscript: `mean_covariance_fit`, `full_gaussian_observed_data_fit`, and `conditional_adequacy_named_variable`.
- Replaced Cholesky whitening as the default global residual basis with numerically safeguarded symmetric/ZCA whitening, $Z=(X-\hat\mu)\hat\Sigma^{-1/2}$.
- Set the eigenvalue threshold to `sqrt(.Machine$double.eps) * max(1, max(abs(eigenvalues)))`; the prototype stops when an implied covariance is not sufficiently positive definite.
- Preserved the former transform as `cholesky_innovations`, labelled as an order-dependent comparison-only representation without original-variable pair labels.
- Added explicit representation labels to global distribution, global dHSIC, marginal distribution, pairwise, and case-distance outputs. Global diagnostics now use `symmetric_zca`; row distances are squared symmetric/ZCA coordinates (equivalently Mahalanobis distances).
- Changed pairwise heatmap axes and pair tables to coordinate-only labels. Named-variable interpretation remains with existing conditional residuals; WP09A did not add conditional tests or plots.
- Revised residual-representation, diagnostic-target, prototype-workflow, simulation caveat, and empirical-example wording to report the failed target/null rather than generic nonlinear SEM misfit.

Acceptance evidence:

- R acceptance checks passed with a positive-definite non-diagonal covariance: symmetric-whitened implied covariance was identity to a maximum absolute error below `1e-10`; a common variable permutation produced the same symmetric/ZCA values under that permutation; and Cholesky innovations changed materially under reordering.
- The API check confirmed three target labels, both residual representations, and symmetric/ZCA row distances.
- `git diff --check` passed.
- `latexmk -pdf -interaction=nonstopmode -halt-on-error misfit.tex` completed successfully. The rebuilt `misfit.pdf` has no fatal errors, unresolved citations/references, or overfull boxes; ordinary underfull-vbox layout warnings remain.

Generated-output consistency:

- `misfit.pdf` was rebuilt and is consistent with the revised text.
- WP09A did not regenerate the simulation tables/figure or the Wage pairwise graphic. Those retained outputs use the earlier representation and the manuscript now labels them as retained illustrative output rather than revised symmetric/ZCA results.

Unresolved issues for later packages:

- WP09B must provide named-variable conditional localization without treating symmetric/ZCA coordinate pairs as original-variable residual pairs.
- A later simulation/output package must regenerate retained simulation and pairwise graphic outputs under the revised representation before they are used as symmetric/ZCA evidence.
- Refit-bootstrap calibration remains unimplemented and untouched.

## Completed: `work_packages/WP09B_conditional_localization.md`

Date: 2026-07-21.

Files changed:

- `sem_misfit_api_prototype.R`
- `test_wp09b_conditional_localization.R`
- `detection_method_notes.md`
- `misfit.tex`
- `misfit.pdf`
- `CURRENT_HANDOFF.md`

WP09B changes:

- Made variable-labelled conditional residuals the default substantive localization output. `conditional_variables` tests $R_j \perp X_{-j}$; `conditional_pairs` reports both directional forms, $R_j \perp X_k$ and $R_k \perp X_j$; optional `conditional_scale` tests $R_j^2 \perp X_{-j}$.
- Each conditional table reports the statistic, distance-correlation effect size when available, raw and Holm-adjusted p-values, conditional-permutation calibration label, exact null, and status. Holm adjustment is applied separately within the variable, directional-pair, and optional scale families.
- Moved symmetric/ZCA pair output to `transformed_pairs` and `sem_misfit_transformed_pairs()`; `sem_misfit_pairs()` now returns directional named-variable localization by default. The transformed heatmap remains explicitly labelled exploratory.
- Added conditional-residual plots against named original predictors, fitted conditional means (`against = "fitted"`), or named moderators (`moderator = ...`).
- Conditional-regression blocks now require positive definiteness under the WP09A relative eigenvalue threshold. The prototype reports the target and conditioning variables in an informative error and applies no hidden regularization.
- Updated the notes and manuscript to state that conditional rejection may reflect mean, scale, or broader conditional-distribution inadequacy, not uniquely a nonlinear mean term.

Acceptance evidence:

- `Rscript test_wp09b_conditional_localization.R` passed. The check uses simulated non-diagonal Gaussian covariance data, confirms direct conditional-regression predictions to below `1e-10`, verifies both directional pair rows, exercises the optional scale table and fitted/moderator plots, and confirms the near-singular conditioning-block error.
- `git diff --check` passed.
- `latexmk -pdf -interaction=nonstopmode -halt-on-error misfit.tex` completed successfully. The rebuilt `misfit.pdf` has no fatal errors, unresolved citations/references, or overfull boxes; ordinary underfull-vbox layout warnings remain.

Generated-output consistency:

- `misfit.pdf` was rebuilt and reflects the conditional-localization contract.
- The retained simulation table/figure and Wage transformed-pair heatmap were not regenerated, by design; they remain labelled as retained transformed-coordinate illustrative output.

Unresolved issues for later packages:

- WP10A must calibrate the conditional variable, directional-pair, and optional scale distance-covariance interfaces with a refit parametric bootstrap. Current p-values are conditional permutation screens and do not account for fitted-parameter uncertainty.
- The conditional API remains limited to the joint Gaussian observed-data interface; it does not yet support case-specific conditional expectations from arbitrary nonlinear SEMs.
- A later output/simulation package must regenerate retained simulation and Wage graphics under the revised contracts before presenting them as current default-localization evidence.

## Calibration Position

## Completed: `work_packages/WP10A_bootstrap_calibration.md`

Date: 2026-07-21.

Files changed:

- `sem_misfit_api_prototype.R`
- `test_wp10a_bootstrap_calibration.R`
- `detection_method_notes.md`
- `CURRENT_HANDOFF.md`

WP10A changes:

- Added distinct `conditional_permutation`, `parametric_fixed`, and `parametric_refit` calibration modes; legacy `permutation` maps to `conditional_permutation`, and legacy `bootstrap` maps to `parametric_refit`.
- `parametric_fixed` simulates Gaussian observed data from the fitted expectation and recomputes diagnostics using that original expectation. `parametric_refit` simulates, refits the identical lavaan model restrictions (or invokes an explicit generic `refit(expectation, data)` callback), and recomputes diagnostics from the refit.
- Parametric modes use raw diagnostic statistics and bootstrap rank p-values, with no nested conditional permutation test in bootstrap replicates. The default remains `B = 999`.
- Variable-, directional-pair-, scale-, and transformed-pair families use each bootstrap replicate's maximum statistic, so their reported bootstrap p-values are family-wise. Global test rows are calibrated against their own retained statistic distribution.
- Added `calibration_details`: requested/successful replicate counts, failure log, failure rule, seed, runtime, rank formula, and retained bootstrap statistics. Failed simulation, refit (including nonconverged or improper lavaan refits), or diagnostic replicates are logged and excluded under the explicit successful-replicate rule.
- Refit support is implemented for lavaan and generic `sem_expectation` objects with supplied `refit` callbacks. OpenMx refitting is not implemented or claimed.

Acceptance evidence:

- `Rscript test_wp10a_bootstrap_calibration.R` passed. It fits a correctly specified Gaussian lavaan model; confirms reproducible fixed/refit results with a fixed seed; verifies the reported rank p-value; confirms the retained conditional-permutation label; exercises the generic callback route; and verifies an intentional refit failure is retained in the failure log while successful replicates are used.
- A default-diagnostic `parametric_refit` lavaan smoke test with `B = 5` also passed, including the global distribution statistic's raw `energy::mvnorm.e` bootstrap path.
- `git diff --check` passed. Its only messages were pre-existing line-ending warnings for unrelated tracked files `SEMmisfit.R` and `SEMmisfit.bib`.

Generated-output consistency:

- WP10A changed no manuscript or simulation-generated outputs. Existing retained simulation and Wage graphics remain unchanged and retain the WP09B consistency caveat.

Unresolved issues for later packages:

- WP10B must use the implemented calibration modes to assess Type I error, runtime, and failure behavior; it owns any simulation/generated-output reruns.
- Parametric bootstrap currently assumes the Gaussian observed-data expectation interface; generic refitting needs a caller-supplied callback, and OpenMx refitting remains unsupported.
- Bootstrap p-values with the default 999 replicates can be costly; no automatic-default policy beyond the explicit user-selected calibration mode has been set.

- Naive asymptotic p-values may be returned only as fast exploratory screens with explicit caveats.
- Permutation tests are useful for residual-dependence diagnostics conditional on fitted residuals, especially pairwise independence tests.
- Parametric bootstrap from the fitted SEM, with refitting on each bootstrap data set, should be the preferred calibration standard for inferential p-values because it accounts for fitted-parameter uncertainty.
- Sample splitting or cross-fitting is a future option if bootstrap calibration is too slow or unstable.

Unresolved calibration concerns:

- Need simulation evidence on Type I error for naive, permutation, and parametric-bootstrap calibration.
- Need runtime/stability evidence for bootstrap-calibrated dHSIC and pairwise tests as the number of observed variables grows.
- Need to decide whether the default package output should compute bootstrap p-values automatically or require an explicit user option because of cost.
- Need to implement a refit-capable bootstrap interface for lavaan and OpenMx before `calibration = "bootstrap"` can be more than an explicit placeholder.

## Manuscript Claims To Keep Softened

- Do not describe HFI or multivariate normality tests as specific detectors of nonlinear misspecification; they detect broad residual distributional deviations.
- State that distributional diagnostics and residual-dependence diagnostics answer different questions.
- Treat current simulation conclusions as provisional until WP05 reruns the simulation using the selected v1 default diagnostics and calibration plan.
- Avoid implying that localization diagnostics identify the correct model revision without follow-up modeling.

## Current Blockers

- No unresolved blocker remains for compiling the current manuscript or aligning it with WP01-WP08 outputs.
- Scientific/manuscript review may still request prose polish, stronger limitations, or a different empirical framing.
- Package API prototype exists, but it is not yet a formal R package and lacks tests.
- Simulation output has been rebuilt for WP05 and integrated, but remains illustrative rather than a full validation study.
- Bootstrap calibration is designed but not implemented.
- The current manuscript/prototype still needs the WP09A--WP11B methodological
  revision before it is treated as submission-ready.

## Files Touched In This Scaffold Or WP01-WP03

- `PROJECT_PLAN.md`
- `CURRENT_HANDOFF.md`
- `detection_method_notes.md`
- `misfit.tex`
- `sem_misfit_api_prototype.R`
- `work_packages/README.md`
- `work_packages/WP01_detection_method.md`
- `work_packages/WP02_software_api.md`
- `work_packages/WP03_manuscript_reframe.md`
- `work_packages/WP04_placeholder_section.md`
- `work_packages/WP05_simulation_rebuild.md`
- `work_packages/WP06_literature_links.md`
- `work_packages/WP07_empirical_example.md`
- `work_packages/WP08_final_integration.md`
- `empirical_example_screen.R`
- `empirical_example_candidate_notes.md`
- `empirical_example_candidates.csv`
- `empirical_example_candidate_tests.csv`
- `empirical_example_candidate_top_pairs.csv`
- `empirical_example_candidate_top_variables.csv`
- `empirical_example_nonlinear_search.R`
- `empirical_example_nonlinear_notes.md`
- `empirical_example_candidate_plots.R`
- `empirical_example_nonlinear_candidates.csv`
- `empirical_example_nonlinear_top_pairs.csv`
- `empirical_example_nonlinear_marginal.csv`
- `empirical_example_nonlinear_smooth_checks.csv`
- `empirical_example_figures/wage_logwage_residual_by_age.png`
- `empirical_example_figures/nhanes_bpsys_residual_by_age.png`
- `empirical_example_figures/nhanes_bpsys_residual_by_bmi.png`
- `empirical_example_wage.R`
- `wage_empirical_summary.csv`
- `wage_empirical_diagnostics.csv`
- `wage_empirical_pairwise.csv`
- `wage_empirical_marginal.csv`
- `wage_empirical_followup.csv`
- `wage_empirical_smooth_checks.csv`
- `empirical_example_figures/wage_linear_residual_by_age.png`
- `empirical_example_figures/wage_spline_residual_by_age.png`

## Exact Next Task

Recommended next task: `work_packages/WP10B_calibration_simulation.md`.

Reason: WP10A now supplies the fixed and refit calibration interfaces, rank
rule, and failure contract. WP10B can evaluate their calibration and runtime
without changing residual representations, localization targets, or the
bootstrap engine.

Handoff rule for WP09A--WP11B: each chat begins by reading this file and its
single named work package, and ends by updating this file with completed
changes, unresolved issues, files touched, whether generated outputs remain
consistent, and the exact next package.
