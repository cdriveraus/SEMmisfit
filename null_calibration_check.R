library(MASS)
library(data.table)
library(goftest)
library(energy)
library(dHSIC)
library(Matrix)

niter <- 500
nperm <- 199
ncores <- max(1, min(10, parallel::detectCores(logical = TRUE) - 2))

combine_pvalues <- function(p_values) {
  p_values <- p_values[is.finite(p_values)]
  if (length(p_values) == 0) return(NA_real_)
  min(1, length(p_values) * min(p_values))
}

make_data <- function(n, nvars) {
  data <- as.data.table(MASS::mvrnorm(
    n = n,
    mu = rep(0.5, nvars),
    Sigma = diag(2, nvars)
  ))
  setnames(data, paste0("V", seq_len(nvars)))
  data[, V3 := 0.3 * V2 + rnorm(n)]
  data
}

fit_observed_moment_model <- function(data) {
  x <- as.matrix(data)
  n <- nrow(x)
  p <- ncol(x)
  mean_hat <- colMeans(x)
  s <- stats::cov(x) * (n - 1) / n
  sigma_hat <- diag(diag(s), p)
  dimnames(sigma_hat) <- dimnames(s)
  allowed_pairs <- list(c(1, 3), c(2, 3))
  for (pair in allowed_pairs) {
    sigma_hat[pair[1], pair[2]] <- s[pair[1], pair[2]]
    sigma_hat[pair[2], pair[1]] <- s[pair[2], pair[1]]
  }
  sigma_hat <- as.matrix(Matrix::nearPD(sigma_hat, keepDiag = TRUE)$mat)
  chi <- n * (
    determinant(sigma_hat, logarithm = TRUE)$modulus[1] +
      sum(diag(s %*% solve(sigma_hat))) -
      determinant(s, logarithm = TRUE)$modulus[1] -
      p
  )
  sem_p <- pchisq(max(chi, 0), df = p * (p - 1) / 2 - length(allowed_pairs), lower.tail = FALSE)
  list(mean = mean_hat, covariance = sigma_hat, sem_p = sem_p)
}

whiten_residuals <- function(data, mean, covariance) {
  raw <- sweep(as.matrix(data), 2, mean, FUN = "-")
  raw %*% solve(chol(covariance))
}

safe_p <- function(expr) tryCatch(expr, error = function(e) NA_real_)

one <- function(n, nvars) {
  data <- make_data(n, nvars)
  fit <- fit_observed_moment_model(data)
  z <- whiten_residuals(data, fit$mean, fit$covariance)
  row_distance <- rowSums(z^2)
  marginal_p <- vapply(seq_len(ncol(z)), function(j) {
    safe_p(goftest::ad.test(z[, j], "pnorm", 0, 1)$p.value)
  }, numeric(1))
  pairwise_dcov_p <- c()
  for (i in seq_len(ncol(z) - 1)) {
    for (j in (i + 1):ncol(z)) {
      pairwise_dcov_p <- c(pairwise_dcov_p, safe_p(energy::dcov.test(z[, i], z[, j], R = nperm)$p.value))
    }
  }
  data.table(
    SEM_p = fit$sem_p,
    Row_distance_AD_p = safe_p(goftest::ad.test(row_distance, "pchisq", df = ncol(z))$p.value),
    Marginal_AD_p = combine_pvalues(marginal_p),
    Energy_MVN_p = safe_p(energy::mvnorm.etest(z, R = nperm)$p.value),
    Global_dHSIC_p = safe_p(dHSIC::dhsic.test(z, matrix.input = TRUE, B = nperm)$p.value),
    Pairwise_dCov_p = combine_pvalues(pairwise_dcov_p)
  )
}

conditions <- CJ(n = c(100, 500), nvars = c(4, 8), iter = seq_len(niter))

cl <- parallel::makeCluster(ncores)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(data.table)
  library(goftest)
  library(energy)
  library(dHSIC)
  library(Matrix)
})
parallel::clusterExport(
  cl,
  c(
    "combine_pvalues",
    "make_data",
    "fit_observed_moment_model",
    "whiten_residuals",
    "safe_p",
    "one",
    "conditions",
    "nperm"
  ),
  envir = environment()
)
parallel::clusterSetRNGStream(cl, 20260710)

message("Running ", nrow(conditions), " no-misfit iterations on ", ncores, " workers")
result_list <- parallel::parLapplyLB(cl, seq_len(nrow(conditions)), function(i) {
  cond <- conditions[i]
  cbind(n = cond$n, nvars = cond$nvars, iter = cond$iter, one(cond$n, cond$nvars))
})
results <- rbindlist(result_list)
save(results, niter, nperm, file = "null_calibration_check_results.rda")

cols <- setdiff(names(results), c("n", "nvars", "iter"))
print(results[, lapply(.SD, function(x) mean(x <= .05, na.rm = TRUE)), by = .(n, nvars), .SDcols = cols])
