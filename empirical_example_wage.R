options(repos = c(CRAN = "https://cloud.r-project.org"))

required <- c("ISLR2", "lavaan", "energy", "dHSIC", "goftest", "ggplot2", "mgcv")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing)
}

source("sem_misfit_api_prototype.R")

set.seed(70713)
B <- 999

data(Wage, package = "ISLR2")
wage <- as.data.frame(Wage)
wage$year_c <- wage$year - min(wage$year)
wage$age_decade <- wage$age / 10
wage$age_c <- wage$age - mean(wage$age)
wage$age_c2 <- wage$age_c^2

# Stage 1: a joint, just-identified SEM for the mean/covariance diagnostic.
core_vars <- c("age_decade", "year_c", "logwage")
wage_core <- wage[, core_vars]
wage_core <- wage_core[stats::complete.cases(wage_core), ]

diagnostic_model <- "
  logwage ~ age_decade + year_c
  age_decade ~~ year_c
"

diagnostic_fit <- lavaan::sem(
  diagnostic_model,
  data = wage_core,
  meanstructure = TRUE,
  estimator = "ML",
  fixed.x = FALSE,
  missing = "listwise"
)

diagnostics <- sem_misfit(
  diagnostic_fit,
  tests = "default",
  calibration = "permutation",
  B = B,
  adjust = "holm",
  seed = 70713
)

tests <- as.data.frame(diagnostics)
pairs <- sem_misfit_pairs(diagnostics)
pairs <- pairs[order(pairs$p_adjusted, -pairs$statistic), ]
variables <- sem_misfit_variables(diagnostics)
variables <- variables[order(variables$p_adjusted, -variables$statistic), ]
cases <- sem_misfit_cases(diagnostics)
fit_measures <- lavaan::fitMeasures(
  diagnostic_fit,
  c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr")
)

# Stage 2: matched conditional SEMs on the same outcome and fixed predictors.
# Keeping age_c2 in both models makes the one-df likelihood-ratio comparison valid.
linear_followup_model <- "
  logwage ~ age_c + 0*age_c2 + year_c
"
quadratic_followup_model <- "
  logwage ~ age_c + age_c2 + year_c
"

linear_sem <- lavaan::sem(
  linear_followup_model,
  data = wage,
  meanstructure = TRUE,
  estimator = "ML",
  fixed.x = TRUE,
  missing = "listwise"
)
quadratic_sem <- lavaan::sem(
  quadratic_followup_model,
  data = wage,
  meanstructure = TRUE,
  estimator = "ML",
  fixed.x = TRUE,
  missing = "listwise"
)

linear_prediction <- as.numeric(lavaan::lavPredictY(linear_sem, newdata = wage)[, "logwage"])
quadratic_prediction <- as.numeric(lavaan::lavPredictY(quadratic_sem, newdata = wage)[, "logwage"])
wage$linear_sem_residual <- wage$logwage - linear_prediction
wage$quadratic_sem_residual <- wage$logwage - quadratic_prediction

sem_lrt <- lavaan::lavTestLRT(quadratic_sem, linear_sem)
quadratic_term <- subset(
  lavaan::parameterEstimates(quadratic_sem),
  lhs == "logwage" & op == "~" & rhs == "age_c2"
)
sem_comparison <- data.frame(
  comparison = "linear versus quadratic age path",
  chisq_difference = sem_lrt[2, "Chisq diff"],
  df_difference = sem_lrt[2, "Df diff"],
  p_value = sem_lrt[2, "Pr(>Chisq)"],
  linear_aic = stats::AIC(linear_sem),
  quadratic_aic = stats::AIC(quadratic_sem),
  delta_aic = stats::AIC(quadratic_sem) - stats::AIC(linear_sem),
  quadratic_estimate = quadratic_term$est,
  quadratic_se = quadratic_term$se,
  quadratic_z = quadratic_term$z,
  quadratic_p = quadratic_term$pvalue,
  stringsAsFactors = FALSE
)

# A spline is retained only as a sensitivity check on the quadratic specification.
gam_age <- mgcv::gam(
  logwage ~ s(age, k = 6) + year_c,
  data = wage,
  method = "REML"
)
wage$spline_residual <- stats::residuals(gam_age)

residual_dependence <- function(residual, label, aic = NA_real_, r_squared = NA_real_) {
  test <- energy::dcov.test(wage$age, residual, R = B)
  data.frame(
    model = label,
    n = length(residual),
    dcor_age_residual = energy::dcor(wage$age, residual),
    dcov_statistic = unname(test$statistic),
    dcov_p = test$p.value,
    residual_sd = stats::sd(residual),
    aic = aic,
    r_squared = r_squared,
    stringsAsFactors = FALSE
  )
}

followup <- rbind(
  residual_dependence(
    wage$linear_sem_residual,
    "linear_conditional_sem",
    stats::AIC(linear_sem),
    unname(lavaan::inspect(linear_sem, "r2")["logwage"])
  ),
  residual_dependence(
    wage$quadratic_sem_residual,
    "quadratic_conditional_sem",
    stats::AIC(quadratic_sem),
    unname(lavaan::inspect(quadratic_sem, "r2")["logwage"])
  ),
  residual_dependence(
    wage$spline_residual,
    "spline_sensitivity",
    stats::AIC(gam_age),
    summary(gam_age)$r.sq
  )
)

gam_summary <- summary(gam_age)
smooth_checks <- data.frame(
  model = "spline_sensitivity",
  smooth = "s(age)",
  edf = gam_summary$s.table[1, "edf"],
  p_value = gam_summary$s.table[1, "p-value"],
  stringsAsFactors = FALSE
)

summary_out <- data.frame(
  dataset = "ISLR2::Wage",
  n = nrow(wage_core),
  variables = paste(core_vars, collapse = ", "),
  model = "logwage ~ age_decade + year_c; age_decade ~~ year_c",
  chisq = unname(fit_measures["chisq"]),
  df = unname(fit_measures["df"]),
  pvalue = unname(fit_measures["pvalue"]),
  cfi = unname(fit_measures["cfi"]),
  rmsea = unname(fit_measures["rmsea"]),
  srmr = unname(fit_measures["srmr"]),
  global_distribution_p = tests$p_value[tests$test == "global_distribution"],
  global_independence_p = tests$p_value[tests$test == "global_independence"],
  pairwise_min_padjusted = tests$p_value[tests$test == "pairwise_independence"],
  marginal_min_padjusted = tests$p_value[tests$test == "marginal_distribution"],
  row_distance_p = tests$p_value[tests$test == "case_distance"],
  stringsAsFactors = FALSE
)

dir.create("empirical_example_figures", showWarnings = FALSE)

pair_plot <- plot(diagnostics, type = "pairwise_heatmap") +
  ggplot2::labs(
    title = "Pairwise localization",
    fill = "dCov"
  ) +
  ggplot2::theme_minimal(base_size = 11)

conditional_plot <- plot(
  diagnostics,
  type = "conditional_residual",
  variable = "logwage",
  against = "age_decade"
) +
  ggplot2::scale_x_continuous(labels = function(x) x * 10) +
  ggplot2::labs(
    x = "Age",
    y = "Conditional residual: log wage",
    title = "Localized nonlinear residual pattern"
  ) +
  ggplot2::theme_minimal(base_size = 11)

quadratic_plot <- ggplot2::ggplot(
  wage,
  ggplot2::aes(.data$age, .data$quadratic_sem_residual)
) +
  ggplot2::geom_point(alpha = 0.18, size = 0.8) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 0.8) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  ggplot2::labs(
    x = "Age",
    y = "Residual log wage",
    title = "Quadratic SEM follow-up"
  ) +
  ggplot2::theme_minimal(base_size = 11)

spline_plot <- ggplot2::ggplot(
  wage,
  ggplot2::aes(.data$age, .data$spline_residual)
) +
  ggplot2::geom_point(alpha = 0.18, size = 0.8) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 0.8) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  ggplot2::labs(
    x = "Age",
    y = "Residual log wage",
    title = "Spline sensitivity check"
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(
  "empirical_example_figures/wage_pairwise_heatmap.png",
  pair_plot,
  width = 6.5,
  height = 4.2,
  dpi = 180
)
ggplot2::ggsave(
  "empirical_example_figures/wage_conditional_residual_by_age.png",
  conditional_plot,
  width = 6.5,
  height = 4.2,
  dpi = 180
)
ggplot2::ggsave(
  "empirical_example_figures/wage_quadratic_residual_by_age.png",
  quadratic_plot,
  width = 6.5,
  height = 4.2,
  dpi = 180
)
ggplot2::ggsave(
  "empirical_example_figures/wage_spline_residual_by_age.png",
  spline_plot,
  width = 6.5,
  height = 4.2,
  dpi = 180
)

write.csv(summary_out, "wage_empirical_summary.csv", row.names = FALSE)
write.csv(tests, "wage_empirical_diagnostics.csv", row.names = FALSE)
write.csv(pairs, "wage_empirical_pairwise.csv", row.names = FALSE)
write.csv(variables, "wage_empirical_marginal.csv", row.names = FALSE)
write.csv(cases, "wage_empirical_cases.csv", row.names = FALSE)
write.csv(followup, "wage_empirical_followup.csv", row.names = FALSE)
write.csv(sem_comparison, "wage_empirical_sem_comparison.csv", row.names = FALSE)
write.csv(smooth_checks, "wage_empirical_smooth_checks.csv", row.names = FALSE)
saveRDS(diagnostics, "wage_empirical_diagnostics.rds")

print(summary_out)
print(tests)
print(pairs)
print(followup)
print(sem_comparison)
print(smooth_checks)
