# WP10B supersedes the earlier illustrative simulation below.  The full
# calibration run lives in WP10B_calibration_simulation.R and is checkpointed.
# ================================================================
# Simulation study for residual diagnostic decomposition in SEM
# ================================================================
#
# This script regenerates simresults.rda, simtable.tex, and sim.pdf.
# The goal is illustrative: compare diagnostics with distinct targets,
# not establish a definitive benchmark of test performance.

required_packages <- c(
  "MASS",
  "data.table",
  "goftest",
  "energy",
  "dHSIC",
  "ggplot2",
  "kableExtra",
  "Matrix"
)

options(repos = c(CRAN = "https://cloud.r-project.org"))

installed_packages <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed_packages) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

set.seed(20260710)

niter <- 500
nperm <- 199
ncores <- max(1, min(20, parallel::detectCores(logical = TRUE) - 2))


# Register parallel backend
cl <- makeCluster(ncores)
registerDoSNOW(cl)

simconditions <- data.table(expand.grid(
  n = c(100, 500),
  nvars = c(4, 8),
  misfit_type = c(
    "none",
    "omitted covariance",
    "omitted interaction",
    "nonlinear curve",
    "skewed marginal",
    "heavy tails",
    "group mean",
    "group variance"
  ),
  stringsAsFactors = FALSE
))

combine_pvalues <- function(p_values, method = "bonferroni") {
  p_values <- p_values[is.finite(p_values)]
  if (length(p_values) == 0) return(NA_real_)
  if (tolower(method) == "fisher") {
    stat <- -2 * sum(log(pmax(p_values, .Machine$double.xmin)))
    return(pchisq(stat, df = 2 * length(p_values), lower.tail = FALSE))
  }
  min(1, length(p_values) * min(p_values))
}

make_data <- function(n, nvars, misfit_type) {
  data <- as.data.table(MASS::mvrnorm(
    n = n,
    mu = rep(0.5, nvars),
    Sigma = diag(2, nvars)
  ))
  setnames(data, paste0("V", seq_len(nvars)))

  v1c <- data$V1 - mean(data$V1)
  v2c <- data$V2 - mean(data$V2)
  eps <- rnorm(n)

  data[, V3 := 0.3 * V2 + eps]

  if (misfit_type == "omitted interaction") {
    data[, V3 := V3 + 0.45 * v1c * v2c]
  }

  if (misfit_type == "nonlinear curve") {
    data[, V3 := V3 + 0.45 * (v1c^2 - mean(v1c^2))]
  }

  if (misfit_type == "skewed marginal") {
    data[, V4 := scale(log1p(exp(V4)))[, 1] * sd(V4) + mean(V4)]
  }

  if (misfit_type == "heavy tails") {
    data[, V4 := 0.5 + sqrt(2) * rt(n, df = 3) / sqrt(3)]
  }

  if (misfit_type == "group mean") {
    idx <- seq_len(floor(n / 2))
    data[idx, c("V2", "V3") := .(V2 + 2, V3 + 2)]
  }

  if (misfit_type == "group variance") {
    idx <- seq_len(floor(n / 2))
    data[idx, c("V2", "V3") := .(0.5 + 2 * (V2 - 0.5), 0.5 + 2 * (V3 - 0.5))]
  }

  data
}

fit_observed_moment_model <- function(data, misfit_type) {
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

  if (misfit_type == "omitted covariance") {
    sigma_hat[2, 3] <- sigma_hat[3, 2] <- 0.05
  }

  sigma_hat <- as.matrix(Matrix::nearPD(sigma_hat, keepDiag = TRUE)$mat)
  chi <- n * (
    determinant(sigma_hat, logarithm = TRUE)$modulus[1] +
      sum(diag(s %*% solve(sigma_hat))) -
      determinant(s, logarithm = TRUE)$modulus[1] -
      p
  )

  n_fixed_cov <- p * (p - 1) / 2 - length(allowed_pairs)
  if (misfit_type == "omitted covariance") {
    n_fixed_cov <- n_fixed_cov + 1
  }
  sem_p <- pchisq(max(chi, 0), df = n_fixed_cov, lower.tail = FALSE)

  list(mean = mean_hat, covariance = sigma_hat, sem_p = sem_p)
}

whiten_residuals <- function(data, mean, covariance) {
  raw <- sweep(as.matrix(data), 2, mean, FUN = "-")
  raw %*% solve(chol(covariance))
}

safe_p <- function(expr) {
  tryCatch(expr, error = function(e) NA_real_)
}

simulate_iteration <- function(cond, iter) {
  data <- make_data(cond$n, cond$nvars, cond$misfit_type)
  fit <- fit_observed_moment_model(data, cond$misfit_type)
  z <- whiten_residuals(data, fit$mean, fit$covariance)
  row_distance <- rowSums(z^2)

  marginal_p <- vapply(seq_len(ncol(z)), function(j) {
    safe_p(goftest::ad.test(z[, j], "pnorm", 0, 1)$p.value)
  }, numeric(1))

  pairwise_dcov_p <- c()
  if (ncol(z) > 1) {
    for (i in seq_len(ncol(z) - 1)) {
      for (j in (i + 1):ncol(z)) {
        pairwise_dcov_p <- c(
          pairwise_dcov_p,
          safe_p(energy::dcov.test(z[, i], z[, j], R = nperm)$p.value)
        )
      }
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

tasks <- simconditions[rep(seq_len(.N), each = niter)]
tasks[, condition := rep(seq_len(nrow(simconditions)), each = niter)]
tasks[, iter := rep(seq_len(niter), times = nrow(simconditions))]

message("Running ", nrow(tasks), " simulation iterations on ", ncores, " workers")

cl <- parallel::makeCluster(ncores)
on.exit(parallel::stopCluster(cl), add = TRUE)
invisible(parallel::clusterEvalQ(cl, {
  library(MASS)
  library(data.table)
  library(goftest)
  library(energy)
  library(dHSIC)
  library(ggplot2)
  library(kableExtra)
  library(Matrix)
}))
parallel::clusterExport(
  cl,
  c(
    "combine_pvalues",
    "make_data",
    "fit_observed_moment_model",
    "whiten_residuals",
    "safe_p",
    "simulate_iteration",
    "tasks",
    "nperm"
  ),
  envir = environment()
)
parallel::clusterSetRNGStream(cl, 20260710)

simresults <- rbindlist(parallel::parLapplyLB(cl, seq_len(nrow(tasks)), function(task_idx) {
  cond <- tasks[task_idx]
  cbind(
    condition = cond$condition,
    n = cond$n,
    nvars = cond$nvars,
    misfit_type = cond$misfit_type,
    simulate_iteration(cond, cond$iter)
  )
}))

save(simresults, simconditions, niter, nperm, file = "simresults.rda")

conditioncols <- c("condition", "n", "nvars", "misfit_type")
outputstatnames <- setdiff(colnames(simresults), conditioncols)

summary_results <- simresults[
  ,
  lapply(.SD, function(x) round(mean(x <= .05, na.rm = TRUE), 3)),
  by = conditioncols,
  .SDcols = outputstatnames
]

print(summary_results)

table <- as.data.frame(summary_results)[, -1]
colnames(table) <- gsub("_p", "", colnames(table), fixed = TRUE)
colnames(table)[colnames(table) == "misfit_type"] <- "Condition"
colnames(table)[colnames(table) == "n"] <- "Subjects"
colnames(table)[colnames(table) == "nvars"] <- "Variables"

sink("simtable.tex")
kableExtra::kable(
  table,
  format = "latex",
  linesep = "",
  booktabs = TRUE,
  digits = 2
) |>
  kableExtra::row_spec(0, bold = TRUE) |>
  kableExtra::kable_styling(latex_options = c("striped"), position = "left")
sink()

summary_long <- data.table::melt(
  summary_results,
  id.vars = c("condition", colnames(simconditions)),
  measure.vars = outputstatnames,
  variable.name = "diagnostic",
  value.name = "detection_rate"
)

summary_long[, diagnostic := factor(
  gsub("_p", "", diagnostic, fixed = TRUE),
  levels = c(
    "SEM",
    "Row_distance_AD",
    "Marginal_AD",
    "Energy_MVN",
    "Global_dHSIC",
    "Pairwise_dCov"
  ),
  labels = c(
    "SEM chi-square",
    "Row-distance AD",
    "Marginal AD",
    "Energy MVN",
    "Global dHSIC",
    "Pairwise dCov"
  )
)]
summary_long[, n := factor(n)]
summary_long[, nvars := factor(nvars)]

p <- ggplot2::ggplot(
  summary_long,
  ggplot2::aes(x = diagnostic, y = detection_rate, fill = diagnostic)
) +
  ggplot2::geom_col(width = 0.78) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(misfit_type),
    cols = ggplot2::vars(interaction(nvars, n, sep = " vars, n=")),
    labeller = ggplot2::label_both
  ) +
  ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    x = NULL,
    y = "Detection rate (p <= .05)",
    fill = "Diagnostic"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    strip.text.y = ggplot2::element_text(angle = 0)
  )

grDevices::pdf("sim.pdf", width = 11, height = 13)
print(p)
grDevices::dev.off()

# ================================================================
# End of simulation study script
# ================================================================
