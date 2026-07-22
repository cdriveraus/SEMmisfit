# Detection Method Notes

Date: 2026-07-09

Scope: WP01 only. These notes set the conceptual diagnostic taxonomy and recommended v1 default tests. They do not change the manuscript, simulation outputs, or package API.

## Diagnostic Taxonomy

The diagnostic workflow should separate the target of a test from the residual representation used by the test. A significant result should be interpreted as evidence against a specific residual null, not as a direct diagnosis of the substantive source of misspecification.

### The Three Diagnostic Targets

The revised prototype uses exactly three labels. **Mean/covariance fit** asks whether the fitted SEM reproduces the observed first and second moments. **Full Gaussian observed-data fit** asks whether observations follow the fitted multivariate Gaussian mean/covariance distribution. **Conditional adequacy for a named variable** asks whether one named variable follows its fitted conditional distribution given named conditioning variables. A test may be a partial probe of one target; its reported null and representation must still be stated explicitly.

### Mean/Covariance Fit

Target: whether the fitted SEM reproduces the observed mean vector and covariance matrix.

Primary residual object: observed minus model-implied mean and covariance structure, already summarized by standard SEM fit tests and residual covariance diagnostics.

Interpretation: detects classical SEM misspecification in first- and second-order moments. It does not directly test nonlinear dependence, distribution shape, heteroscedasticity, or subgroup structure unless those features change the mean or covariance enough to affect the fitted SEM.

Recommended role: keep as the baseline comparator, not as the new nonlinear diagnostic.

### Full Gaussian Observed-Data Fit: Marginal Distribution Probe

Target: whether each whitened residual coordinate follows the expected univariate standard normal distribution under the fitted model.

Primary residual object: columns of whitened observed residuals.

Interpretation: sensitive to skewness, kurtosis, tail behavior, outliers, and variable-specific nonnormality. This is a distributional diagnostic. It should not be described as specifically detecting nonlinear structural misspecification because marginal nonnormality can arise without residual dependence.

Recommended role: optional or secondary v1 diagnostic, useful for flagging distributional violations and localizing them by variable.

### Full Gaussian Observed-Data Fit: Joint Distribution Probe

Target: whether the full whitened residual vector follows the expected multivariate normal distribution with zero mean and identity covariance.

Primary residual object: matrix of whitened observed residuals.

Interpretation: sensitive to broad departures from the fitted Gaussian observed-data distribution, including marginal nonnormality, tail dependence, mixture structure, outliers, and unmodeled dependence. This is a global distributional diagnostic, not a pure nonlinear-association diagnostic.

Recommended role: v1 default as an omnibus residual distribution check, reported separately from residual-dependence tests.

### Full Gaussian Observed-Data Fit: Residual-Dependence Probe

Target: whether residual coordinates are mutually independent after whitening by the fitted model-implied covariance.

Primary residual object: columns of whitened observed residuals.

Interpretation: if the fitted SEM fully captures the observed mean and covariance and the residual vector is Gaussian, whitening should remove linear dependence. Remaining dependence among residual coordinates suggests unmodeled higher-order association, nonlinear dependence, mixture/subgroup structure, or other joint structure not represented by the model. These tests are less directly sensitive to purely marginal nonnormality than distributional tests, but they are not source-specific.

Recommended role: central v1 diagnostic target. Use a joint independence test as the main global nonlinear-dependence diagnostic and pairwise independence tests for localization.

### Heteroscedasticity and Subgroup Structure

Target: whether residual scale or residual distance varies systematically across observations, variables, fitted values, observed covariates, or latent subgroups.

Primary residual objects: rowwise squared residual distance, conditional residuals, residual absolute values or squares, fitted values/predictions, moderators, or observed grouping variables.

Interpretation: can detect omitted interactions, mixture structure, subgroup mean/variance differences, and outlier clusters. However, a scalar fourth-moment index can also react to nonnormality and outliers, so it is not a clean heteroscedasticity test unless calibrated and localized against a specified moderator or fitted-value axis.

Recommended role: do not include the current HFI-style scalar index as a default p-value. Include rowwise distance and moderator plots/summaries as exploratory localization, and simulate formal heteroscedasticity tests before promoting them.

### Localization by Variable, Pair, Case, or Moderator

Target: where a global diagnostic signal appears.

Primary residual objects:

- variable-level marginal residuals,
- pairwise residual dependence tests,
- case-level Mahalanobis/rowwise squared residual distances,
- residual scale by fitted values or moderators,
- subgroup comparisons when a candidate grouping variable is known.

Interpretation: localization is diagnostic triage. It identifies variables, pairs, cases, or moderators worth inspecting, but it should not be treated as confirmatory model selection without multiplicity control and follow-up modeling.

Recommended role: expose alongside global diagnostics. V1 should return local statistics and adjusted p-values where possible, but frame them as localization rather than automatic modification advice.

## Residual Representations

### Symmetric/ZCA Whitened Observed Residuals

Definition: for each row, subtract the model-implied observed mean and multiply by the symmetric inverse square root of the model-implied observed covariance, $Z=(X-\hat\mu)\hat\Sigma^{-1/2}$. The prototype uses an eigendecomposition and rejects covariance matrices whose eigenvalues do not exceed $\sqrt{\epsilon_{\mathrm{machine}}}\max(1,\max|\lambda|)$.

Expected behavior under the fitted Gaussian SEM: approximately independent standard normal coordinates, subject to fitted-parameter uncertainty.

Strengths:

- generic interface for SEM software that can provide implied observed means and covariances;
- directly supports marginal normality, multivariate normality, rowwise distance, and residual-dependence tests;
- clear null under the fitted observed-data model.

Limitations:

- estimated means/covariances are fit to the same data, so naive p-values can be anti-conservative or otherwise miscalibrated;
- symmetric/ZCA coordinates are equivariant to a common reordering of the data and implied covariance, but they remain transformed coordinates rather than original-variable residuals;
- coordinate-level marginal and pairwise outputs must not be read as variable-labelled localization.

Recommended v1 role: default global representation for the full Gaussian observed-data target, including global distribution, dHSIC, marginal-distribution, and case-distance diagnostics.

### Cholesky Innovations (Optional Comparison Only)

The prototype preserves the existing upper-Cholesky transform as `cholesky_innovations`. It is an order-dependent innovation representation, not the default global basis. Its columns have no original-variable pair labels and must not be used to localize a named-variable signal. Squared row distances agree with symmetric/ZCA distances despite this coordinate-order sensitivity.

### Conditional Residuals for Each Variable Given the Rest

Definition: for each variable, subtract its fitted conditional expectation given all other observed variables under the model-implied joint Gaussian distribution, then standardize by the fitted conditional variance.

Expected behavior under the fitted Gaussian SEM: each conditional residual is standard normal and independent of the conditioning set in population.

Strengths:

- easier variable-level interpretation than arbitrary whitened coordinates;
- useful for plots against fitted conditional predictions, moderators, and other variables;
- natural representation for localizing which observed variable has unexplained structure.

Limitations:

- conditional residuals for different target variables are not mutually independent objects in the same way as whitened coordinates;
- running global independence tests across all conditional residual columns has a less transparent null;
- calibration still depends on fitted parameters.

Recommended v1 role: named-variable localization and plotting, not the main global test basis. The prototype tests $R_j \perp X_{-j}$ with multivariate distance covariance and directional $R_j \perp X_k$ tests with the same statistic. `conditional_permutation` remains an explicitly conditional exploratory screen. WP10A additionally implements `parametric_fixed` (simulate from the fitted expectation and reuse it) and `parametric_refit` (simulate, refit the same restrictions, then recompute); the latter is implemented for lavaan and generic expectations with an explicit refit callback. An optional conditional-scale screen tests $R_j^2 \perp X_{-j}$. Bootstrap p-values use the retained-statistic rank and bootstrap maxima within each selected localization family for family-wise inference.

Numerical handling: each conditioning covariance block must be positive definite under the same relative eigenvalue threshold used for symmetric whitening. The prototype stops with the target and conditioning-variable names when that condition fails; it deliberately applies no hidden regularization.

### Rowwise Squared Residual Distances

Definition: sum of squared symmetric/ZCA residual coordinates, equivalent to a Mahalanobis distance from the fitted mean/covariance.

Expected behavior under a known Gaussian mean/covariance: chi-square with degrees of freedom equal to the number of observed variables.

Strengths:

- simple case-level misfit measure;
- useful for outlier detection, mixture/subgroup hints, and distributional omnibus checks;
- invariant to rotations of the whitened residual vector when computed as a Mahalanobis distance.

Limitations:

- collapses direction and variable/pair information;
- deviations from chi-square can reflect marginal nonnormality, tail behavior, subgroup structure, heteroscedasticity, or dependence;
- chi-square p-values ignore that the mean/covariance are estimated.

Recommended v1 role: default descriptive/localization output; include a calibrated distributional test only with bootstrap calibration or clear caveat.

### Pairwise Residuals

Definition: the default substantive pairs are directional pairs of a variable-specific conditional residual and an original labelled predictor or moderator. Symmetric/ZCA coordinate pairs are retained only as explicitly transformed exploratory output.

Expected behavior under the fitted Gaussian SEM: for each named target, $R_j$ is independent of its conditioning set and therefore of each named predictor $X_k$.

Strengths:

- preserves the direction and labels of $R_j$ versus $X_k$ for follow-up;
- supports distance covariance, pairwise dHSIC/HSIC, HHG-style tests, and scatter/residual plots;
- closer to actionable model revision than a single global rejection.

Limitations:

- many comparisons as the number of variables grows;
- pairwise independence does not guarantee joint independence;
- $R_j \perp X_k$ is directional in its representation even when both directional tests are reported;
- multiple variable, pair, and optional scale screens require separate multiplicity control.

Recommended v1 role: default substantive localization, with Holm adjustment separately for variable-level, directional-pair, and optional scale families. A rejection can reflect conditional mean, conditional scale, or another conditional-distribution failure; it does not uniquely identify a nonlinear mean term.

## Calibration Plan

### Naive Asymptotic P-Values

Use only as fast screening statistics or clearly labeled exploratory output. They generally assume a fixed null distribution and do not account for SEM parameters estimated from the same data. For residual diagnostics based on fitted means/covariances, this can distort Type I error.

Naive p-values are most defensible for preliminary ranking or visual triage, not final claims about calibrated rejection rates.

### Permutation Tests

Use for residual-dependence tests where the null is exchangeability/independence among residual coordinates. For pairwise tests, permuting one residual coordinate relative to another gives a clear conditional null of pairwise independence. For joint dHSIC-style tests, column-wise or variable-block permutations can approximate the mutual-independence null.

Permutation alone does not account for model refitting uncertainty. It conditions on the fitted residuals and therefore tests whether the residual matrix contains dependence beyond its fitted marginal distributions. This is useful for dependence diagnostics but should not be presented as a fully unconditional SEM goodness-of-fit p-value.

### Parametric Bootstrap from the Fitted SEM

Recommended calibration standard for v1 when computation allows.

Procedure:

1. Fit the SEM to the observed data.
2. Compute all diagnostics from the observed fitted residuals.
3. Simulate bootstrap data from the fitted SEM-implied observed distribution, using the fitted mean and covariance for the initial continuous single-group complete-data scope.
4. Refit the same SEM to each bootstrap data set.
5. Recompute the same diagnostics after refitting.
6. Compare observed diagnostics with the refit bootstrap distribution.

This directly addresses fitted-parameter uncertainty and the effect of model estimation on residual distributions. It also gives a common calibration layer across distributional and dependence diagnostics.

Implementation contract: `parametric_fixed` and `parametric_refit` use $B=999$ by default and do not run nested permutation tests in bootstrap replicates. The prototype records the requested and successful replicate counts, seed, elapsed runtime, and a failure log. Failed simulation, refit, or diagnostic replicates are retained in that log and excluded only under the stated successful-replicate rule. `parametric_refit` currently supports lavaan and generic expectation objects only when the caller supplies a `refit(expectation, data)` callback; OpenMx refitting remains unsupported.

### Sample Splitting or Cross-Fitting

Useful if bootstrap cost is too high or if overfitting/calibration problems persist. Fit the SEM on one part of the data, compute residual diagnostics on held-out rows, and repeat across folds for stability.

Benefits:

- diagnostics are computed on data not used to estimate the implied moments;
- naive null approximations become closer to valid.

Costs:

- lower effective sample size per diagnostic;
- more implementation complexity;
- less standard for SEM users.

Recommended v1 role: not the default, but keep as a future calibration option and a simulation condition in WP05 if bootstrap calibration is too slow.

## Recommended V1 Default Diagnostic Set

The package default should return both global tests and localization outputs. Defaults should be narrow, interpretable, and calibrated by parametric bootstrap when the user asks for inferential p-values.

### Baseline SEM Moment Fit

Statistic: standard SEM chi-square test and/or available global fit indices from the fitted SEM software.

Null hypothesis: the model-implied mean/covariance structure is sufficient for the observed first- and second-order moments under the estimator's assumptions.

Interpretation: classical moment/covariance fit baseline. It should be reported separately from nonlinear residual diagnostics.

V1 status: default comparator, not a new residual diagnostic.

### Global Residual Distribution Test

Statistic: energy-distance multivariate goodness-of-fit test or equivalent omnibus test comparing symmetric/ZCA residual vectors to the standard multivariate normal reference; alternatively, a rowwise Mahalanobis-distance Anderson-Darling test can be provided as a simpler scalar distribution check.

Null hypothesis: the full Gaussian observed-data target holds: fitted symmetric/ZCA residual vectors are standard multivariate normal, after accounting for fitting through the chosen calibration method.

Interpretation: broad distributional misfit. A rejection can reflect marginal nonnormality, tail behavior, mixtures, outliers, heteroscedasticity, or residual dependence. It is not specific to nonlinear dependence.

V1 status: default global distributional diagnostic. Prefer bootstrap-calibrated p-values; naive p-values only with a caveat.

### Global Residual Independence Test

Statistic: joint dHSIC on symmetric/ZCA residual coordinates, with permutation or parametric-bootstrap calibration.

Null hypothesis: the dHSIC probe finds no residual-coordinate dependence within the full Gaussian observed-data target.

Interpretation: evidence of unmodeled joint residual dependence, including nonlinear association, interactions, mixtures, or other higher-order structure not captured by the fitted mean/covariance. This diagnostic is conceptually distinct from marginal and joint distributional tests.

V1 status: primary nonlinear-dependence diagnostic.

### Pairwise Residual Independence Localization

Statistic: directional distance covariance between $R_j$ and original variable $X_k$, with multiplicity-adjusted p-values and distance-correlation effect sizes. Symmetric/ZCA coordinate-pair statistics are available only as transformed exploratory output.

Null hypothesis: $R_j \perp X_k$ within the named-variable conditional-adequacy target.

Interpretation: a rejection can reflect conditional mean, scale, or broader conditional-distribution inadequacy, rather than uniquely an omitted nonlinear mean term. Use original-variable conditional-residual plots against predictors, fitted values, or prespecified moderators.

V1 status: default substantive localization output when the number of observed variables is modest; for high-dimensional models, return top-ranked pairs and adjusted p-values.

### Marginal Residual Normality Localization

Statistic: univariate Anderson-Darling tests of each symmetric/ZCA residual coordinate against standard normality, with multiplicity adjustment.

Null hypothesis: each symmetric/ZCA residual coordinate follows a standard normal distribution within the full Gaussian observed-data target.

Interpretation: variable/coordinate-level distributional violation, not residual dependence. Useful for identifying skewness, heavy tails, and outlying coordinates.

V1 status: secondary default or optional default. Include if clearly labeled as distributional localization.

### Case-Level Residual Distance

Statistic: rowwise squared Mahalanobis distance from the fitted observed mean/covariance, plus chi-square reference quantiles or bootstrap reference quantiles.

Null hypothesis for inferential use: the case-level distances follow the fitted-model reference distribution after accounting for parameter fitting.

Interpretation: case-level extremeness or subgroup/mixure hints. It does not identify which dependency is missing by itself.

V1 status: default descriptive output; inferential p-values only when bootstrap-calibrated.

## Diagnostics to Exclude From V1 Defaults

### Current HFI Scalar as a P-Value

Reason to exclude: the current script returns a transformed fourth-moment-style index with an ad hoc offset and threshold. It is sensitive to heteroscedasticity, kurtosis, outliers, and nonnormal errors, but the implemented value is not a calibrated p-value. It should not be mixed with p-values in summary tables.

Possible future role: descriptive heteroscedasticity/tail diagnostic after its null, scaling, threshold, and calibration are re-established.

### Mutual Information With Discretization

Reason to exclude: discretization choices affect the statistic, power, and Type I error. The current aggregate variants also have less transparent null hypotheses than dHSIC or distance covariance.

Possible future role: experimental comparison in simulations if a stable binning rule and calibration scheme are specified.

### Aggregate Residual-Rest Mutual Information

Reason to exclude: summing the remaining residual columns before testing mutual information with a target coordinate creates a hard-to-interpret composite. Rejections are difficult to map back to a substantive residual-dependence pattern.

Possible future role: not needed for v1 if joint dHSIC plus pairwise localization are available.

### Levene-Style Tests Using Residual-Derived Groups

Reason to exclude: grouping by residuals or arbitrary cuts can create circularity and unclear nulls. Heteroscedasticity tests are useful only when the moderator, fitted value, group, or conditioning axis is specified before testing.

Possible future role: moderator-specific residual-scale diagnostics.

### Overlap Coefficients Against Simulated References

Reason to exclude: less standard inferential target and less direct interpretability than energy-distance, Anderson-Darling, or bootstrap reference distributions.

Possible future role: descriptive distribution similarity measure if users need effect-size-like summaries.

### Naive Chi-Square Test of Rowwise Distances as a Confirmatory P-Value

Reason to exclude: the chi-square reference assumes known mean/covariance. The SEM-implied moments are estimated, and the model can absorb some sample deviations. Use row distances descriptively or calibrate by parametric bootstrap.

## Implications for WP02

WP02 should expose a small set of package defaults:

- `sem_fit`: baseline SEM fit statistic from the fitted model object when available.
- `resid_global_distribution`: bootstrap-calibrated global distributional diagnostic for whitened residuals.
- `resid_global_independence`: joint dHSIC-style global residual-independence diagnostic.
- `resid_pairwise_independence`: pairwise distance covariance or HSIC/dHSIC localization with multiplicity adjustment.
- `resid_marginal_distribution`: optional/default-secondary univariate residual normality localization.
- `case_distance`: rowwise squared Mahalanobis distances with reference quantiles.
- `conditional_residuals`: variable-labeled localization residuals for plots and moderator checks.

The API should distinguish statistic values, naive p-values, permutation p-values, and parametric-bootstrap p-values. It should not label HFI-like index values as p-values.

## Implications for WP05

WP05 should simulate at least these questions:

- Type I error of naive, permutation, and parametric-bootstrap calibration under correctly specified Gaussian SEMs.
- Type I error under marginal nonnormality that does not introduce residual dependence.
- Power for classical covariance misfit, omitted interaction/nonlinear dependence, subgroup mean shifts, subgroup variance shifts, and outliers/heavy tails.
- Whether global dHSIC remains relatively insensitive to pure marginal nonnormality compared with global distributional tests.
- Whether pairwise distance covariance or pairwise HSIC/dHSIC better localizes omitted interactions.
- Whether rowwise distance diagnostics identify subgroup and outlier structure without being overinterpreted as nonlinear association.
- Runtime and stability as sample size and number of observed variables increase.

## Manuscript Claims to Soften

The manuscript should avoid saying that HFI or multivariate normality tests specifically detect nonlinear misspecification. They detect broad residual distributional deviations. The clean distinction should be:

- distributional diagnostics test whether whitened residuals match the full fitted Gaussian residual distribution;
- residual-dependence diagnostics test whether residual coordinates remain independent after whitening;
- localization diagnostics suggest where to inspect, but do not identify the correct model revision on their own.

Claims based on the current simulation should remain provisional until WP05 reruns the simulation with the v1 diagnostic set and bootstrap/permutation calibration.
