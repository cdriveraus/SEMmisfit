source("sem_misfit_api_prototype.R")
stopifnot(requireNamespace("lavaan", quietly = TRUE), requireNamespace("energy", quietly = TRUE))

set.seed(20260721)
n <- 120L
dat <- data.frame(x = rnorm(n), y = rnorm(n), z = rnorm(n))
dat$y <- .45 * dat$x + dat$y
dat$z <- -.30 * dat$x + .25 * dat$y + dat$z
model <- "y ~ x\nz ~ x + y\nx ~~ x"
fit <- lavaan::sem(model, data = dat, meanstructure = TRUE)

run_mode <- function(mode) sem_misfit(
  fit, tests = c("conditional_variables", "conditional_pairs", "conditional_scale"),
  calibration = mode, B = 19L, seed = 4242
)

fixed_a <- run_mode("parametric_fixed")
fixed_b <- run_mode("parametric_fixed")
refit_a <- run_mode("parametric_refit")
refit_b <- run_mode("parametric_refit")
for (x in list(fixed_a, refit_a)) {
  stopifnot(x$calibration_details$B_successful > 0L)
  stopifnot(all(x$conditional_variables$calibration == x$calibration))
  stopifnot(all(x$conditional_pairs$status == "bootstrap_familywise_rank"))
  boot <- x$calibration_details$bootstrap_statistics$conditional_variables
  expected <- bootstrap_rank_p(x$conditional_variables$statistic[1], boot)
  stopifnot(isTRUE(all.equal(x$conditional_variables$p_value[1], expected)))
}
stopifnot(identical(fixed_a$conditional_pairs$p_value, fixed_b$conditional_pairs$p_value))
stopifnot(identical(refit_a$conditional_pairs$p_value, refit_b$conditional_pairs$p_value))

perm <- sem_misfit(fit, tests = "conditional_variables", calibration = "conditional_permutation", B = 19L, seed = 4242)
stopifnot(all(perm$conditional_variables$calibration == "conditional_permutation"))

generic <- as_sem_expectation(list(
  data = dat, mean = colMeans(dat), covariance = cov(dat), variables = names(dat),
  simulate = function(expectation, n) as.data.frame(MASS::mvrnorm(n, expectation$mean, expectation$covariance)),
  refit = function(expectation, data) as_sem_expectation(expectation, data = data)
))
generic_refit <- sem_misfit(generic, tests = "conditional_variables", calibration = "parametric_refit", B = 9L, seed = 91)
stopifnot(generic_refit$calibration_details$B_successful == 9L)

failure_counter <- 0L
generic_with_failure <- generic
generic_with_failure$refit <- function(expectation, data) {
  failure_counter <<- failure_counter + 1L
  if (failure_counter == 1L) stop("intentional test refit failure")
  as_sem_expectation(expectation, data = data)
}
failure_logged <- sem_misfit(generic_with_failure, tests = "conditional_variables", calibration = "parametric_refit", B = 5L, seed = 92)
stopifnot(failure_logged$calibration_details$B_successful == 4L)
stopifnot(nrow(failure_logged$calibration_details$failures) == 1L)
stopifnot(grepl("intentional test refit failure", failure_logged$calibration_details$failures$reason[1]))

cat("WP10A bootstrap calibration checks passed.\n")
