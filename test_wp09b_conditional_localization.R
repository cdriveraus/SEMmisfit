source("sem_misfit_api_prototype.R")

set.seed(90902)
variables <- c("outcome", "predictor_a", "predictor_b")
mu <- c(outcome = 1.5, predictor_a = -0.5, predictor_b = 2)
sigma <- matrix(
  c(2.0, 0.6, -0.3,
    0.6, 1.5, 0.4,
   -0.3, 0.4, 1.2),
  nrow = 3, dimnames = list(variables, variables)
)
data <- as.data.frame(MASS::mvrnorm(250, mu = mu, Sigma = sigma))
names(data) <- variables

residuals <- compute_sem_residuals(
  as_sem_expectation(list(data = data, mean = mu, covariance = sigma, variables = variables))
)
rest <- c("predictor_a", "predictor_b")
direct_prediction <- mu["outcome"] +
  sigma["outcome", rest, drop = FALSE] %*% solve(sigma[rest, rest, drop = FALSE]) %*%
  t(sweep(as.matrix(data[rest]), 2, mu[rest], "-"))
stopifnot(max(abs(residuals$conditional$prediction[, "outcome"] - as.vector(direct_prediction))) < 1e-10)

diagnostics <- sem_misfit(
  list(data = data, mean = mu, covariance = sigma, variables = variables),
  B = 49, seed = 90903
)
stopifnot(
  all(c("conditional_variables", "conditional_pairs", "conditional_scale") %in% names(diagnostics)),
  nrow(diagnostics$conditional_variables) == 3L,
  nrow(diagnostics$conditional_pairs) == 6L,
  nrow(diagnostics$conditional_scale) == 0L,
  all(c("statistic", "effect_size", "p_value", "calibration", "null_hypothesis") %in% names(diagnostics$conditional_variables)),
  all(c("residual_variable", "predictor") %in% names(diagnostics$conditional_pairs)),
  any(diagnostics$conditional_pairs$residual_variable == "outcome" & diagnostics$conditional_pairs$predictor == "predictor_a"),
  any(diagnostics$conditional_pairs$residual_variable == "predictor_a" & diagnostics$conditional_pairs$predictor == "outcome"),
  !identical(sem_misfit_pairs(diagnostics), sem_misfit_transformed_pairs(diagnostics))
)

scale_diagnostics <- sem_misfit(
  list(data = data, mean = mu, covariance = sigma, variables = variables),
  tests = "conditional_scale", B = 49, seed = 90904
)
stopifnot(nrow(scale_diagnostics$conditional_scale) == 3L)

fitted_plot <- plot(diagnostics, type = "conditional_residual", variable = "outcome", against = "fitted")
moderator_plot <- plot(diagnostics, type = "conditional_residual", variable = "outcome", moderator = "predictor_a")
stopifnot(inherits(fitted_plot, "ggplot"), inherits(moderator_plot, "ggplot"))

near_singular <- sigma
near_singular["predictor_b", ] <- near_singular["predictor_a", ]
near_singular[, "predictor_b"] <- near_singular[, "predictor_a"]
near_singular["predictor_b", "predictor_b"] <- near_singular["predictor_a", "predictor_a"]
near_singular_error <- tryCatch({
  compute_conditional_residuals(as.matrix(data), mu, near_singular, variables)
  NULL
}, error = conditionMessage)
stopifnot(!is.null(near_singular_error), grepl("Conditioning covariance block|conditional variance", near_singular_error))

cat("WP09B conditional-localization checks passed.\n")
