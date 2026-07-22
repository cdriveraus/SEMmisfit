# Literature Notes

Date: 2026-07-10

Scope: WP06 literature links only. These notes map concrete citation support for the current manuscript and identify remaining references to read before finalizing claims.

## Existing Bibliography Entries To Use More Strongly

- `bollen1989structural`: classical SEM foundation and the mean/covariance target of SEM fit.
- `bentler1990comparative`, `hu1999cutoff`: global fit-index context; use only for ordinary SEM fit, not residual-distribution diagnostics.
- `micceri1989unicorn`: empirical motivation for nonnormal psychological data.
- `gerhard2017fit`: HFI as the directly related nonlinear-SEM residual-variance proposal; cite when explaining why HFI is relevant but not a default calibrated p-value here.
- `mcdonald1962general`: early nonlinear factor-analysis precedent; useful for positioning nonlinear follow-up modeling as established rather than invented here.
- `rizzo2016energy`: energy distance for multivariate distribution checks and distance-covariance family methods.
- `pfister2018kernelbased`: dHSIC, permutation/bootstrap/gamma calibration, and joint residual-independence tests.
- `karch2024pearsons`: psychological-methods review of modern nonparametric independence tests, including distance correlation and HHG-style tests.
- `brandmaier2013structural`: SEM trees for subgroup and covariate-driven heterogeneity follow-up.
- `briley2015nonparametric`: LOSEM / moderator-varying SEM follow-up for smooth covariate patterns.
- `muthen2012bayesian`: Bayesian SEM, approximate zero restrictions, and posterior predictive checking.

## Recent SEM-Specific Layer Added After WP06 Review

- `irmer2024estimating`: recent nonlinear SEM / moderation power and software context; supports the claim that nonlinear SEM remains an active methodological area.
- `cox2025comparing`: recent comparison of Bayesian and structural-after-measurement approaches for latent-interaction SEM with complex data structures.
- `jia2026obtaining`: recent robust practical fit-index work for multiply imputed nonnormal SEM data; supports the distinction between robust SEM fit handling and residual-diagnostic decomposition.
- `garnierVillarreal2020adapting`: Bayesian SEM fit-index adaptation; supports the Bayesian fit/posterior-predictive follow-up sentence.
- `arnold2021score`: score-guided SEM trees as a recent development beyond the original SEM tree citation.
- `silvaDiaz2025evaluation`: SEM forests for omitted influential covariates / specification search; supports adding forests alongside SEM trees as heterogeneity follow-up tools.

## References Added In WP06

### Classical SEM fit, residual diagnostics, and calibration

- `bollen1992bootstrapping`: supports the claim that SEM fit and residual diagnostics based on fitted models need bootstrap-style calibration when inferential p-values are desired.
- `satorra1994corrections`: supports robust corrections for SEM test statistics and standard errors under nonnormality. This should be described as robust inference for moment-structured SEM, not as a residual-dependence diagnostic.

### HFI and nonlinear SEM

- No new HFI reference was needed because `gerhard2017fit` was already present and is now cited directly in the HFI paragraph.
- `klein2000maximum`: added for latent moderated structural equations / latent interaction follow-up when residual-dependence diagnostics suggest higher-order association.

### Mixture and subgroup follow-up

- `lubke2005investigating`: added for factor mixture models as a follow-up when rowwise distances, residual-scale patterns, or heatmaps suggest discrete heterogeneity.

### Recent nonlinear, fit, and heterogeneity SEM literature

- `irmer2024estimating`, `cox2025comparing`, `jia2026obtaining`, `garnierVillarreal2020adapting`, `arnold2021score`, and `silvaDiaz2025evaluation` were added after a follow-up check that the SEM-specific literature layer was too historically anchored.
- These references are used in the manuscript introduction, SEM fit paragraph, and follow-up-model paragraph. They should not be expanded into a broad review unless the manuscript scope changes.

## Missing References Still Worth Reading Before Final Claims

### Robust/non-normal SEM and heavy tails

- Need a focused read on robust ML / sandwich / scaled-test options beyond `satorra1994corrections` and `jia2026obtaining`, especially if the manuscript makes more detailed claims about heavy-tailed SEM estimation. Candidate area: Yuan and Bentler robust mean/covariance structure analysis and more recent lavaan-oriented robust test-statistic comparisons.
- Need a citation specifically for heavy-tailed or t-distribution SEM if the follow-up paragraph is expanded beyond generic robust corrections.

### Residual whitening, Mahalanobis distances, and case diagnostics

- The current manuscript uses basic multivariate-normal and Mahalanobis-distance logic, but it still lacks a direct citation for Mahalanobis distance / multivariate outlier diagnostics. Add one only if the row-distance section becomes more than a brief derivation.
- Case-distance p-values should remain caveated until refit-bootstrap calibration is implemented or simulated.

### dHSIC/HSIC, distance covariance, energy tests, HHG

- `rizzo2016energy`, `pfister2018kernelbased`, and `karch2024pearsons` are enough for the current manuscript scope.
- Add original distance covariance / HSIC / HHG papers only if the methods section gives technical definitions or contrasts among those tests.

### SEM trees, LOSEM, and model-based recursive partitioning

- `brandmaier2013structural`, `briley2015nonparametric`, `arnold2021score`, and `silvaDiaz2025evaluation` cover the current follow-up claims.
- A general model-based recursive partitioning citation would be useful only if the manuscript broadens beyond SEM trees and forests to the general statistical framework.

### Bayesian posterior predictive checks

- `muthen2012bayesian` and `garnierVillarreal2020adapting` support the current BSEM/Bayesian-fit/posterior-predictive sentence.
- Add a broader posterior predictive checking citation only if posterior predictive checking becomes a central calibration recommendation.

### Empirical example source

- No empirical example source has been selected. WP07 should identify a data source that can illustrate at least one residual-distribution signal and one residual-dependence/localization signal without overclaiming model repair.

## Manuscript Locations Patched In WP06

- Introduction: positioned novelty as integration and workflow rather than invention of all diagnostics.
- Introduction: added recent SEM-specific citations for nonlinear SEM, latent interactions, robust fit indices, Bayesian fit indices, SEM trees, and SEM forests.
- Introduction follow-up paragraph: added citations for nonlinear factor analysis, latent interactions, mixture/latent-class follow-up, SEM trees, LOSEM, and Bayesian SEM.
- SEM & Model Fit: added robust-correction/bootstrap citations while preserving the point that classical fit indices target mean/covariance reproduction.
- Rowwise residual distances: added bootstrap citation for inferential calibration.
- Distance/kernel independence paragraph: added energy/dHSIC/independence-test citations and named HHG-style tests.
- Heteroscedasticity paragraph: cited HFI directly and kept the calibration caveat.
- Simulation interpretation: cited Bollen-Stine and dHSIC calibration sources when calling for refit bootstrap or stronger calibration.
- Follow-up-model paragraph: tied robust corrections, recent robust fit-index work, latent interactions, nonlinear factor analysis, mixture SEM, SEM trees/forests, LOSEM, and Bayesian SEM fit indices to specific diagnostic signals.

## Remaining Concrete Citation Gaps In `misfit.tex`

- Add a direct heavy-tailed/t-distribution SEM reference if the paper recommends heavy-tailed SEM as more than a possible follow-up.
- Add an original Mahalanobis/outlier-diagnostics reference if rowwise distances become a major method rather than a descriptive residual summary.
- Add original distance covariance / HSIC / HHG references only if the manuscript explains these methods technically.
- Add empirical-example citations in WP07 once the example data set is chosen.
