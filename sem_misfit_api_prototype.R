# Prototype API for residual-based SEM misfit diagnostics.
#
# Scope: continuous observed variables, single group, complete/listwise data.
# This file is intentionally package-like but not yet a full R package.

as_sem_expectation <- function(object, data = NULL, ...) {
  UseMethod("as_sem_expectation")
}

as_sem_expectation.sem_expectation <- function(object, data = NULL, ...) {
  if (!is.null(data)) {
    object$data <- as.data.frame(data)
    object$n <- nrow(object$data)
  }
  validate_sem_expectation(object)
}

as_sem_expectation.list <- function(object, data = NULL, ...) {
  if (!is.null(data)) {
    object$data <- data
  }

  aliases <- list(
    mean = c("mean", "means", "mu", "implied_mean", "expMean"),
    covariance = c("covariance", "cov", "sigma", "implied_cov", "expCov"),
    variables = c("variables", "vars", "var_names", "names", "manifestVars"),
    sem_fit = c("sem_fit", "fit", "fit_measures")
  )

  pick <- function(x, names) {
    hit <- intersect(names, names(x))
    if (length(hit)) x[[hit[1]]] else NULL
  }

  out <- list(
    data = pick(object, "data"),
    mean = pick(object, aliases$mean),
    covariance = pick(object, aliases$covariance),
    variables = pick(object, aliases$variables),
    n = object$n %||% NA_integer_,
    group = object$group %||% NULL,
    source = object$source %||% "list",
    sem_fit = pick(object, aliases$sem_fit),
    simulate = object$simulate %||% NULL,
    refit = object$refit %||% NULL
  )
  class(out) <- "sem_expectation"
  validate_sem_expectation(out)
}

as_sem_expectation.lavaan <- function(object, data = NULL, ...) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("The lavaan package is required to read lavaan objects.", call. = FALSE)
  }

  observed <- data %||% lavaan::lavInspect(object, "data")
  implied <- lavaan::lavInspect(object, "implied")

  if (is.list(observed) && !is.data.frame(observed)) {
    if (length(observed) != 1L) {
      stop("v1 supports single-group lavaan fits only.", call. = FALSE)
    }
    observed <- observed[[1L]]
  }
  if (is.list(implied) && !is.null(implied[[1L]]) && is.list(implied[[1L]])) {
    if (length(implied) != 1L) {
      stop("v1 supports single-group lavaan fits only.", call. = FALSE)
    }
    implied <- implied[[1L]]
  }

  implied_cov <- implied$cov %||% implied$cov.ov
  implied_mean <- implied$mean %||% implied$mean.ov
  if (is.null(implied_mean)) {
    implied_mean <- rep(0, ncol(implied_cov))
  }

  out <- list(
    data = as.data.frame(observed),
    mean = as.numeric(implied_mean),
    covariance = as.matrix(implied_cov),
    variables = colnames(implied_cov) %||% names(implied_mean),
    n = nrow(observed),
    group = NULL,
    source = "lavaan",
    fit = object,
    sem_fit = tryCatch(
      lavaan::fitMeasures(object, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr")),
      error = function(e) NULL
    )
  )
  class(out) <- "sem_expectation"
  validate_sem_expectation(out)
}

as_sem_expectation.MxModel <- function(object, data = NULL, ...) {
  attrs <- find_openmx_fitfunction_expectations(object)
  observed <- data %||% openmx_observed_data(object)
  if (is.null(observed)) {
    stop("OpenMx data could not be found; pass data explicitly.", call. = FALSE)
  }

  out <- list(
    data = as.data.frame(observed),
    mean = as.numeric(attrs$expMean),
    covariance = as.matrix(attrs$expCov),
    variables = colnames(attrs$expCov) %||% colnames(observed),
    n = nrow(observed),
    group = NULL,
    source = "OpenMx",
    sem_fit = openmx_fit_summary(object)
  )
  class(out) <- "sem_expectation"
  validate_sem_expectation(out)
}

sem_misfit <- function(object,
                       data = NULL,
                       tests = "default",
                       calibration = c("auto", "conditional_permutation", "parametric_fixed", "parametric_refit", "naive", "permutation", "bootstrap", "none"),
                       B = 999,
                       adjust = "holm",
                       seed = NULL,
                       ...) {
  calibration <- match.arg(calibration)
  expectation <- if (inherits(object, "sem_expectation")) {
    as_sem_expectation(object, data = data)
  } else {
    as_sem_expectation(object, data = data, ...)
  }

  tests <- expand_sem_misfit_tests(tests)
  residuals <- compute_sem_residuals(expectation)
  diagnostics <- compute_sem_diagnostics(
    expectation = expectation,
    residuals = residuals,
    tests = tests,
    calibration = calibration,
    B = B,
    adjust = adjust,
    seed = seed
  )
  if (calibration %in% c("parametric_fixed", "parametric_refit", "bootstrap")) {
    diagnostics <- bootstrap_calibrate_diagnostics(
      expectation, tests, calibration = calibration, B = B, adjust = adjust,
      seed = seed, observed = diagnostics
    )
  }

  out <- list(
    call = match.call(),
    expectation = expectation,
    scope = c(
      "continuous observed variables",
      "single group",
      "complete/listwise data",
      "model-implied observed mean and covariance"
    ),
    calibration = diagnostics$calibration,
    targets = sem_misfit_targets(),
    residuals = residuals,
    tests = diagnostics$tests,
    # `pairs` remains available only as transformed-coordinate exploration.
    transformed_pairs = diagnostics$transformed_pairs,
    pairs = diagnostics$conditional_pairs,
    variables = diagnostics$variables,
    conditional_variables = diagnostics$conditional_variables,
    conditional_pairs = diagnostics$conditional_pairs,
    conditional_scale = diagnostics$conditional_scale,
    cases = diagnostics$cases,
    plot_data = diagnostics$plot_data,
    calibration_details = diagnostics$calibration_details %||% NULL
  )
  class(out) <- "sem_misfit"
  out
}

print.sem_misfit <- function(x, ...) {
  cat("SEM residual misfit diagnostics\n")
  cat("Source:", x$expectation$source, "\n")
  cat("Scope:", paste(x$scope, collapse = "; "), "\n")
  cat("Calibration:", x$calibration, "\n\n")
  print(x$tests[, c("test", "target", "representation", "statistic", "p_value", "method", "status")], row.names = FALSE)
  invisible(x)
}

as.data.frame.sem_misfit <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$tests
}

sem_misfit_pairs <- function(x) {
  stopifnot(inherits(x, "sem_misfit"))
  # Default substantive localization is directional and variable-labelled.
  x$pairs
}

sem_misfit_transformed_pairs <- function(x) {
  stopifnot(inherits(x, "sem_misfit"))
  x$transformed_pairs
}

sem_misfit_variables <- function(x) {
  stopifnot(inherits(x, "sem_misfit"))
  x$variables
}

sem_misfit_conditional_variables <- function(x) {
  stopifnot(inherits(x, "sem_misfit"))
  x$conditional_variables
}

sem_misfit_conditional_scale <- function(x) {
  stopifnot(inherits(x, "sem_misfit"))
  x$conditional_scale
}

sem_misfit_cases <- function(x) {
  stopifnot(inherits(x, "sem_misfit"))
  x$cases
}

plot.sem_misfit <- function(x,
                            type = c("pairwise_heatmap", "rowwise_misfit", "conditional_residual"),
                            variable = NULL,
                            against = NULL,
                            moderator = NULL,
                            ...) {
  type <- match.arg(type)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for plot methods.", call. = FALSE)
  }

  if (type == "pairwise_heatmap") {
    pairs <- x$transformed_pairs
    if (is.null(pairs) || !nrow(pairs)) {
      stop("No pairwise diagnostic data are available.", call. = FALSE)
    }
    mirror <- pairs
    mirror$var1 <- pairs$var2
    mirror$var2 <- pairs$var1
    pairs <- rbind(pairs, mirror)
    coordinate_levels <- colnames(x$residuals$symmetric_zca)
    pairs$var1 <- factor(pairs$var1, levels = coordinate_levels)
    pairs$var2 <- factor(pairs$var2, levels = coordinate_levels)
    pairs$label <- paste0(
      "dCov = ", formatC(pairs$statistic, digits = 2, format = "f"),
      "\n", "adjusted p = ", format.pval(pairs$p_adjusted, digits = 2)
    )
    return(
      ggplot2::ggplot(pairs, ggplot2::aes(.data$var1, .data$var2, fill = .data$statistic)) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 3) +
        ggplot2::coord_equal() +
        ggplot2::scale_x_discrete(drop = FALSE) +
        ggplot2::scale_y_discrete(drop = FALSE) +
        ggplot2::labs(
          x = NULL, y = NULL, fill = "Statistic",
          subtitle = "Symmetric-ZCA coordinates; not original-variable residual pairs"
        )
    )
  }

  if (type == "conditional_residual") {
    available <- x$expectation$variables
    if (is.null(variable) || !variable %in% available) {
      stop(
        "'variable' must name one modeled variable: ",
        paste(available, collapse = ", "),
        call. = FALSE
      )
    }
    if (!is.null(moderator)) {
      if (!moderator %in% available) {
        stop("'moderator' must name one modeled variable: ", paste(available, collapse = ", "), call. = FALSE)
      }
      x_values <- x$expectation$data[[moderator]]
      x_label <- moderator
    } else if (identical(against, "fitted")) {
      x_values <- x$residuals$conditional$prediction[, variable]
      x_label <- paste("Fitted conditional mean:", variable)
    } else if (is.null(against) || !against %in% available) {
      stop(
        "'against' must name a modeled predictor or be 'fitted'; alternatively supply 'moderator'. ",
        "Modeled variables: ", paste(available, collapse = ", "),
        call. = FALSE
      )
    } else if (identical(variable, against)) {
      stop("'variable' and 'against' must name different variables.", call. = FALSE)
    } else {
      x_values <- x$expectation$data[[against]]
      x_label <- against
    }

    plot_data <- data.frame(
      against = x_values,
      residual = x$residuals$conditional$standardized[, variable],
      stringsAsFactors = FALSE
    )
    return(
      ggplot2::ggplot(plot_data, ggplot2::aes(.data$against, .data$residual)) +
        ggplot2::geom_point(alpha = 0.18, size = 0.8) +
        ggplot2::geom_smooth(
          method = "loess",
          formula = y ~ x,
          se = TRUE,
          linewidth = 0.8
        ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
        ggplot2::labs(
          x = x_label,
          y = paste("Standardized conditional residual:", variable)
        )
    )
  }

  cases <- x$cases
  ggplot2::ggplot(cases, ggplot2::aes(.data$case, .data$row_distance)) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Case", y = "Squared symmetric-ZCA residual distance")
}

validate_sem_expectation <- function(x) {
  required <- c("data", "mean", "covariance")
  missing <- required[vapply(x[required], is.null, logical(1))]
  if (length(missing)) {
    stop("Missing expectation fields: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  x$data <- as.data.frame(x$data)
  x$covariance <- as.matrix(x$covariance)
  x$mean <- as.numeric(x$mean)

  if (!all(vapply(x$data, is.numeric, logical(1)))) {
    stop("v1 supports continuous numeric observed variables only.", call. = FALSE)
  }
  if (anyNA(x$data)) {
    x$data <- stats::na.omit(x$data)
    x$n <- nrow(x$data)
  }
  if (length(x$mean) != ncol(x$covariance)) {
    stop("Implied mean length must match implied covariance dimension.", call. = FALSE)
  }
  if (nrow(x$covariance) != ncol(x$covariance)) {
    stop("Implied covariance must be square.", call. = FALSE)
  }
  if (ncol(x$data) != length(x$mean)) {
    stop("Observed data columns must match implied mean/covariance dimension.", call. = FALSE)
  }

  if (is.null(x$variables) || length(x$variables) != length(x$mean)) {
    x$variables <- colnames(x$data) %||% paste0("V", seq_along(x$mean))
  }
  colnames(x$data) <- x$variables
  names(x$mean) <- x$variables
  dimnames(x$covariance) <- list(x$variables, x$variables)
  x$n <- nrow(x$data)
  x$group <- x$group %||% NULL
  x$source <- x$source %||% "unknown"
  class(x) <- "sem_expectation"
  x
}

compute_sem_residuals <- function(expectation) {
  data <- as.matrix(expectation$data)
  mu <- expectation$mean
  sigma <- expectation$covariance
  variables <- expectation$variables

  raw <- sweep(data, 2, mu, FUN = "-")
  chol_upper <- chol(sigma)
  # Base R's chol() returns upper R with sigma = t(R) %*% R.
  # This comparison representation is order-dependent and must not be used
  # for variable-labelled localization.
  cholesky_innovations <- t(backsolve(chol_upper, t(raw), transpose = TRUE))

  symmetric_zca <- symmetric_zca_whiten(raw, sigma)
  colnames(raw) <- variables
  colnames(symmetric_zca) <- paste0("zca_coordinate_", seq_len(ncol(symmetric_zca)))
  colnames(cholesky_innovations) <- paste0("cholesky_innovation_", seq_len(ncol(cholesky_innovations)))

  conditional <- compute_conditional_residuals(data, mu, sigma, variables)
  row_distance <- rowSums(symmetric_zca^2)

  list(
    raw = raw,
    symmetric_zca = symmetric_zca,
    cholesky_innovations = cholesky_innovations,
    # Backward-compatible alias. All global diagnostics use symmetric_zca.
    whitened = symmetric_zca,
    conditional = conditional,
    row_distance = row_distance,
    representations = list(
      global = "symmetric_zca",
      cholesky_innovations = "optional_comparison_only_order_dependent",
      eigenvalue_tolerance = sqrt(.Machine$double.eps)
    )
  )
}

symmetric_zca_whiten <- function(raw, sigma,
                                 eigenvalue_tolerance = sqrt(.Machine$double.eps)) {
  sigma_symmetric <- (sigma + t(sigma)) / 2
  decomposition <- eigen(sigma_symmetric, symmetric = TRUE)
  scale <- max(1, max(abs(decomposition$values)))
  minimum <- eigenvalue_tolerance * scale
  if (any(!is.finite(decomposition$values)) || any(decomposition$values <= minimum)) {
    stop(
      "Model-implied covariance is not sufficiently positive definite for symmetric/ZCA whitening ",
      "(minimum eigenvalue must exceed ", format(minimum, scientific = TRUE), ").",
      call. = FALSE
    )
  }
  inverse_square_root <- decomposition$vectors %*%
    diag(1 / sqrt(decomposition$values), nrow = length(decomposition$values)) %*%
    t(decomposition$vectors)
  raw %*% inverse_square_root
}

compute_conditional_residuals <- function(data, mu, sigma, variables) {
  n <- nrow(data)
  p <- ncol(data)
  prediction <- residual <- standardized <- matrix(NA_real_, n, p)
  colnames(prediction) <- colnames(residual) <- colnames(standardized) <- variables

  for (j in seq_len(p)) {
    rest <- setdiff(seq_len(p), j)
    sigma_jj <- sigma[j, j, drop = TRUE]
    if (!length(rest)) {
      if (!is.finite(sigma_jj) || sigma_jj <= 0) {
        stop("Conditional residual for '", variables[j], "' requires a positive implied variance.", call. = FALSE)
      }
      prediction[, j] <- mu[j]
      residual[, j] <- data[, j] - prediction[, j]
      standardized[, j] <- residual[, j] / sqrt(sigma_jj)
      next
    }
    sigma_j_rest <- sigma[j, rest, drop = FALSE]
    sigma_rest_rest <- sigma[rest, rest, drop = FALSE]
    sigma_inv <- conditional_block_inverse(sigma_rest_rest, variables[rest], variables[j])
    cond_var <- as.numeric(sigma_jj - sigma_j_rest %*% sigma_inv %*% t(sigma_j_rest))
    if (!is.finite(cond_var) || cond_var <= 0) {
      stop("Implied conditional variance for '", variables[j], "' is not positive after conditioning on ",
           paste(variables[rest], collapse = ", "), ".", call. = FALSE)
    }
    prediction[, j] <- mu[j] + as.vector(sigma_j_rest %*% sigma_inv %*% t(data[, rest, drop = FALSE] - rep(mu[rest], each = n)))
    residual[, j] <- data[, j] - prediction[, j]
    standardized[, j] <- residual[, j] / sqrt(cond_var)
  }

  list(
    prediction = prediction,
    residual = residual,
    standardized = standardized
  )
}

conditional_block_inverse <- function(block, conditioning_variables, target_variable,
                                      eigenvalue_tolerance = sqrt(.Machine$double.eps)) {
  block <- (block + t(block)) / 2
  values <- eigen(block, symmetric = TRUE, only.values = TRUE)$values
  threshold <- eigenvalue_tolerance * max(1, max(abs(values)))
  if (any(!is.finite(values)) || any(values <= threshold)) {
    stop(
      "Conditioning covariance block for '", target_variable, "' is not sufficiently positive definite ",
      "(conditioning variables: ", paste(conditioning_variables, collapse = ", "), "; minimum eigenvalue must exceed ",
      format(threshold, scientific = TRUE), "). No regularization is applied in v1.",
      call. = FALSE
    )
  }
  solve(block)
}

compute_sem_diagnostics <- function(expectation,
                                    residuals,
                                    tests,
                                    calibration,
                                    B,
                                    adjust,
                                    seed) {
  if (!is.null(seed)) set.seed(seed)
  effective_calibration <- resolve_calibration(calibration)
  test_rows <- list()
  transformed_pair_rows <- data.frame()
  variable_rows <- data.frame()
  conditional_variable_rows <- conditional_table_template("variable")
  conditional_pair_rows <- conditional_table_template("pair")
  conditional_scale_rows <- conditional_table_template("scale")

  if ("sem_fit" %in% tests) {
    test_rows$sem_fit <- sem_fit_row(expectation)
  }
  if ("case_distance" %in% tests) {
    test_rows$case_distance <- row_distance_test_row(residuals$row_distance, ncol(residuals$symmetric_zca))
  }
  if ("marginal_distribution" %in% tests) {
    marginal <- marginal_distribution_rows(residuals$symmetric_zca, adjust = adjust)
    test_rows$marginal_distribution <- marginal$global
    variable_rows <- rbind(variable_rows, marginal$variables)
  }
  if ("global_distribution" %in% tests) {
    test_rows$global_distribution <- global_distribution_row(residuals$symmetric_zca, B = B,
                                                              calibration = effective_calibration)
  }
  if ("global_independence" %in% tests) {
    test_rows$global_independence <- global_independence_row(
      residuals$symmetric_zca,
      calibration = effective_calibration,
      B = B
    )
  }
  if ("pairwise_independence" %in% tests) {
    transformed_pair_rows <- pairwise_independence_rows(
      residuals$symmetric_zca,
      calibration = effective_calibration,
      B = B,
      adjust = adjust
    )
    test_rows$pairwise_independence <- pairwise_global_row(transformed_pair_rows)
  }
  if ("conditional_variables" %in% tests) {
    conditional_variable_rows <- conditional_variable_test_rows(
      residuals$conditional$standardized, expectation$data, calibration = effective_calibration,
      B = B, adjust = adjust
    )
  }
  if ("conditional_pairs" %in% tests) {
    conditional_pair_rows <- conditional_pair_test_rows(
      residuals$conditional$standardized, expectation$data, calibration = effective_calibration,
      B = B, adjust = adjust
    )
  }
  if ("conditional_scale" %in% tests) {
    conditional_scale_rows <- conditional_scale_test_rows(
      residuals$conditional$standardized, expectation$data, calibration = effective_calibration,
      B = B, adjust = adjust
    )
  }

  tests_df <- do.call(rbind, test_rows)
  rownames(tests_df) <- NULL

  cases <- data.frame(
    case = seq_along(residuals$row_distance),
    row_distance = residuals$row_distance,
    reference_quantile = stats::qchisq(
      stats::ppoints(length(residuals$row_distance)),
      df = ncol(residuals$symmetric_zca)
    )
  )

  list(
    calibration = effective_calibration,
    tests = tests_df,
    transformed_pairs = transformed_pair_rows,
    variables = variable_rows,
    conditional_variables = conditional_variable_rows,
    conditional_pairs = conditional_pair_rows,
    conditional_scale = conditional_scale_rows,
    cases = cases,
    plot_data = list(
      pairwise_heatmap = transformed_pair_rows,
      conditional_pairs = conditional_pair_rows,
      rowwise_misfit = cases
    )
  )
}

expand_sem_misfit_tests <- function(tests) {
  if (identical(tests, "default")) {
    return(c(
      "sem_fit",
      "global_distribution",
      "global_independence",
      "pairwise_independence",
      "conditional_variables",
      "conditional_pairs",
      "marginal_distribution",
      "case_distance"
    ))
  }
  match.arg(
    tests,
    choices = c(
      "sem_fit",
      "global_distribution",
      "global_independence",
      "pairwise_independence",
      "conditional_variables",
      "conditional_pairs",
      "conditional_scale",
      "marginal_distribution",
      "case_distance"
    ),
    several.ok = TRUE
  )
}

resolve_calibration <- function(calibration) {
  if (calibration == "auto") {
    return("permutation_for_dependence_naive_for_distribution")
  }
  if (calibration == "bootstrap") {
    return("parametric_refit")
  }
  if (calibration == "permutation") return("conditional_permutation")
  calibration
}

# Parametric-bootstrap calibration deliberately reuses *statistics*, not the
# permutation p-values used for the conditional exploratory screens.  For a
# table family the bootstrap maximum is used, yielding family-wise p-values.
bootstrap_calibrate_diagnostics <- function(expectation, tests, calibration, B,
                                            adjust, seed, observed) {
  mode <- resolve_calibration(calibration)
  if (!mode %in% c("parametric_fixed", "parametric_refit")) stop("Unknown bootstrap mode.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  started <- proc.time()[["elapsed"]]
  replicate_diagnostics <- vector("list", B)
  failures <- vector("list", 0L)

  for (b in seq_len(B)) {
    simulated <- tryCatch(simulate_sem_expectation(expectation), error = function(e) e)
    if (inherits(simulated, "error")) {
      failures[[length(failures) + 1L]] <- bootstrap_failure_row(b, "simulate", simulated)
      next
    }
    boot_expectation <- tryCatch({
      if (mode == "parametric_fixed") as_sem_expectation(expectation, data = simulated) else refit_sem_expectation(expectation, simulated)
    }, error = function(e) e)
    if (inherits(boot_expectation, "error")) {
      failures[[length(failures) + 1L]] <- bootstrap_failure_row(b, "refit", boot_expectation)
      next
    }
    replicate_diagnostics[[b]] <- tryCatch(
      compute_sem_diagnostics(boot_expectation, compute_sem_residuals(boot_expectation), tests,
                              calibration = mode, B = 0L, adjust = adjust, seed = NULL),
      error = function(e) e
    )
    if (inherits(replicate_diagnostics[[b]], "error")) {
      failures[[length(failures) + 1L]] <- bootstrap_failure_row(b, "diagnostic", replicate_diagnostics[[b]])
      replicate_diagnostics[[b]] <- NULL
    }
  }
  retained <- Filter(Negate(is.null), replicate_diagnostics)
  if (!length(retained)) stop("No successful parametric-bootstrap replicates; inspect calibration_details$failures.", call. = FALSE)
  observed <- bootstrap_apply_p_values(observed, retained, mode)
  observed$calibration <- mode
  observed$calibration_details <- list(
    mode = mode, B_requested = B, B_successful = length(retained),
    failures = if (length(failures)) do.call(rbind, failures) else data.frame(),
    failure_rule = "Simulation, refit, or diagnostic failures are logged and excluded; inference is reported only from successful replicates.",
    seed = seed, runtime_seconds = unname(proc.time()[["elapsed"]] - started),
    p_value = "(1 + number of retained bootstrap statistics >= observed statistic) / (1 + successful replicates)",
    bootstrap_statistics = bootstrap_statistic_store(retained)
  )
  observed
}

bootstrap_statistic_store <- function(replicates) {
  out <- list(tests = lapply(replicates, function(x) stats::setNames(x$tests$statistic, x$tests$test)))
  for (nm in c("transformed_pairs", "conditional_variables", "conditional_pairs", "conditional_scale")) {
    out[[nm]] <- vapply(replicates, function(x) max_or_na(x[[nm]]$statistic), numeric(1))
  }
  out
}

simulate_sem_expectation <- function(expectation) {
  if (is.function(expectation$simulate)) return(as.data.frame(expectation$simulate(expectation, expectation$n)))
  z <- matrix(stats::rnorm(expectation$n * length(expectation$mean)), ncol = length(expectation$mean))
  out <- sweep(z %*% chol(expectation$covariance), 2L, expectation$mean, "+")
  colnames(out) <- expectation$variables
  as.data.frame(out)
}

refit_sem_expectation <- function(expectation, data) {
  if (is.function(expectation$refit)) {
    fitted <- expectation$refit(expectation, data)
    return(if (inherits(fitted, "sem_expectation")) as_sem_expectation(fitted) else as_sem_expectation(fitted, data = data))
  }
  if (identical(expectation$source, "lavaan") && !is.null(expectation$fit)) {
    fitted <- lavaan::update(expectation$fit, data = data, do.fit = TRUE)
    if (!lavaan::lavInspect(fitted, "converged")) stop("lavaan refit did not converge.", call. = FALSE)
    if (!isTRUE(lavaan::lavInspect(fitted, "post.check"))) {
      stop("lavaan refit is improper (failed post-estimation admissibility check).", call. = FALSE)
    }
    return(as_sem_expectation(fitted))
  }
  stop("parametric_refit requires a lavaan fit or a generic sem_expectation$refit callback.", call. = FALSE)
}

bootstrap_failure_row <- function(replicate, stage, error) data.frame(
  replicate = replicate, stage = stage, reason = conditionMessage(error), stringsAsFactors = FALSE
)

bootstrap_rank_p <- function(observed_statistic, bootstrap_statistics) {
  bootstrap_statistics <- bootstrap_statistics[is.finite(bootstrap_statistics)]
  if (!is.finite(observed_statistic) || !length(bootstrap_statistics)) return(NA_real_)
  (1 + sum(bootstrap_statistics >= observed_statistic)) / (1 + length(bootstrap_statistics))
}

bootstrap_apply_p_values <- function(observed, replicates, mode) {
  for (test in observed$tests$test) {
    j <- match(test, observed$tests$test)
    boot <- vapply(replicates, function(x) x$tests$statistic[match(test, x$tests$test)], numeric(1))
    observed$tests$p_value[j] <- bootstrap_rank_p(observed$tests$statistic[j], boot)
    observed$tests$method[j] <- paste0(observed$tests$method[j], "; parametric-bootstrap rank")
    observed$tests$status[j] <- "bootstrap_rank"
  }
  for (nm in c("transformed_pairs", "conditional_variables", "conditional_pairs", "conditional_scale")) {
    tab <- observed[[nm]]
    if (is.null(tab) || !nrow(tab)) next
    boot_max <- vapply(replicates, function(x) max_or_na(x[[nm]]$statistic), numeric(1))
    tab$p_value <- vapply(tab$statistic, bootstrap_rank_p, numeric(1), bootstrap_statistics = boot_max)
    tab$p_adjusted <- tab$p_value
    tab$calibration <- mode
    tab$status <- "bootstrap_familywise_rank"
    tab$method <- paste0(tab$method, "; parametric-bootstrap family maximum")
    observed[[nm]] <- tab
  }
  observed
}

sem_fit_row <- function(expectation) {
  fit <- expectation$sem_fit
  p_value <- NA_real_
  statistic <- NA_real_
  if (!is.null(fit)) {
    statistic <- unname(fit[intersect(names(fit), c("chisq", "Chi", "chi"))[1]]) %||% NA_real_
    p_value <- unname(fit[intersect(names(fit), c("pvalue", "p", "P"))[1]]) %||% NA_real_
  }
  diagnostic_row(
    test = "sem_fit",
    target = "mean_covariance_fit",
    statistic = statistic,
    p_value = p_value,
    method = "source_fit_measure",
    status = if (is.null(fit)) "unavailable" else "ok",
    representation = "observed_mean_covariance"
  )
}

row_distance_test_row <- function(row_distance, df) {
  if (requireNamespace("goftest", quietly = TRUE)) {
    test <- goftest::ad.test(row_distance, "pchisq", df = df)
    statistic <- unname(test$statistic)
    p_value <- test$p.value
    method <- "Anderson-Darling chi-square reference"
    status <- "naive_p_value"
  } else {
    statistic <- mean(row_distance)
    p_value <- NA_real_
    method <- "mean row distance; install goftest for AD p-value"
    status <- "statistic_only"
  }
  diagnostic_row(
    "case_distance", "full_gaussian_observed_data_fit", statistic, p_value,
    method, status, representation = "symmetric_zca"
  )
}

marginal_distribution_rows <- function(whitened, adjust) {
  variables <- colnames(whitened)
  rows <- lapply(seq_along(variables), function(j) {
    if (requireNamespace("goftest", quietly = TRUE)) {
      test <- goftest::ad.test(whitened[, j], "pnorm", 0, 1)
      statistic <- unname(test$statistic)
      p_value <- test$p.value
      method <- "Anderson-Darling standard normal reference"
      status <- "naive_p_value"
    } else {
      statistic <- mean(whitened[, j]^3)
      p_value <- NA_real_
      method <- "sample skewness proxy; install goftest for AD p-value"
      status <- "statistic_only"
    }
    data.frame(
      variable = variables[j],
      representation = "symmetric_zca_coordinate",
      statistic = statistic,
      p_value = p_value,
      method = method,
      status = status
    )
  })
  variable_df <- do.call(rbind, rows)
  variable_df$p_adjusted <- stats::p.adjust(variable_df$p_value, method = adjust)

  diagnostic <- diagnostic_row(
    test = "marginal_distribution",
    target = "full_gaussian_observed_data_fit",
    statistic = min_or_na(variable_df$statistic),
    p_value = min_or_na(variable_df$p_adjusted),
    method = paste("minimum adjusted variable-level p-value;", adjust),
    status = unique(variable_df$status)[1],
    representation = "symmetric_zca"
  )
  list(global = diagnostic, variables = variable_df)
}

global_distribution_row <- function(whitened, B, calibration = "auto") {
  if (requireNamespace("energy", quietly = TRUE)) {
    if (calibration %in% c("parametric_fixed", "parametric_refit")) {
      return(diagnostic_row(
        test = "global_distribution",
        target = "full_gaussian_observed_data_fit",
        statistic = unname(energy::mvnorm.e(whitened)), p_value = NA_real_,
        method = "energy::mvnorm.e (parametric-bootstrap rank)",
        status = "bootstrap_rank_pending", representation = "symmetric_zca"
      ))
    }
    test <- energy::mvnorm.etest(whitened, R = B)
    return(diagnostic_row(
      test = "global_distribution",
      target = "full_gaussian_observed_data_fit",
      statistic = unname(test$statistic),
      p_value = test$p.value,
      method = "energy::mvnorm.etest against multivariate normality",
      status = "naive_or_internal_resampling",
      representation = "symmetric_zca"
    ))
  }

  diagnostic_row(
    test = "global_distribution",
    target = "full_gaussian_observed_data_fit",
    statistic = mean(rowSums(whitened^2)),
    p_value = NA_real_,
    method = "mean row distance; install energy for mvnorm.etest",
    status = "statistic_only",
    representation = "symmetric_zca"
  )
}

global_independence_row <- function(whitened, calibration, B) {
  if (requireNamespace("dHSIC", quietly = TRUE)) {
    method <- if (grepl("permutation|auto", calibration)) "permutation" else "gamma"
    test <- dHSIC::dhsic.test(
      whitened,
      matrix.input = TRUE,
      B = if (method == "permutation") B else NULL,
      method = method
    )
    return(diagnostic_row(
      test = "global_independence",
      target = "full_gaussian_observed_data_fit",
      statistic = unname(test$statistic),
      p_value = test$p.value,
      method = paste("dHSIC", method),
      status = "ok",
      representation = "symmetric_zca"
    ))
  }

  diagnostic_row(
    test = "global_independence",
    target = "full_gaussian_observed_data_fit",
    statistic = NA_real_,
    p_value = NA_real_,
    method = "install dHSIC for joint independence test",
    status = "unavailable",
    representation = "symmetric_zca"
  )
}

conditional_table_template <- function(kind) {
  if (kind == "variable") {
    return(data.frame(
      variable = character(), conditioning_variables = character(), representation = character(),
      statistic = numeric(), effect_size = numeric(), p_value = numeric(), p_adjusted = numeric(),
      method = character(), calibration = character(), null_hypothesis = character(), status = character(),
      stringsAsFactors = FALSE
    ))
  }
  if (kind == "pair") {
    return(data.frame(
      residual_variable = character(), predictor = character(), representation = character(),
      statistic = numeric(), effect_size = numeric(), p_value = numeric(), p_adjusted = numeric(),
      method = character(), calibration = character(), null_hypothesis = character(), status = character(),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    variable = character(), conditioning_variables = character(), representation = character(),
    statistic = numeric(), effect_size = numeric(), p_value = numeric(), p_adjusted = numeric(),
    method = character(), calibration = character(), null_hypothesis = character(), status = character(),
    stringsAsFactors = FALSE
  )
}

conditional_dependence_test <- function(response, conditioning, calibration, B) {
  conditioning <- as.matrix(conditioning)
  if (requireNamespace("energy", quietly = TRUE)) {
    if (calibration %in% c("parametric_fixed", "parametric_refit")) {
      statistic <- unname(energy::dcov(response, conditioning))
      effect_size <- tryCatch(as.numeric(energy::dcor(response, conditioning)), error = function(e) NA_real_)
      return(list(statistic = statistic, effect_size = effect_size, p_value = NA_real_,
                  method = "energy::dcov (parametric-bootstrap rank)", calibration = calibration,
                  status = "bootstrap_rank_pending"))
    }
    test <- energy::dcov.test(response, conditioning, R = B)
    effect_size <- tryCatch(
      as.numeric(energy::dcor(response, conditioning)),
      error = function(e) NA_real_
    )
    return(list(
      statistic = unname(test$statistic), effect_size = effect_size, p_value = test$p.value,
      method = "energy::dcov.test", calibration = if (grepl("permutation|auto", calibration)) "conditional_permutation" else "permutation_used",
      status = "ok"
    ))
  }

  correlations <- vapply(seq_len(ncol(conditioning)), function(k) {
    abs(stats::cor(response, conditioning[, k], method = "spearman"))
  }, numeric(1))
  list(
    statistic = max_or_na(correlations), effect_size = max_or_na(correlations), p_value = NA_real_,
    method = "maximum Spearman correlation proxy; install energy for distance covariance",
    calibration = "statistic_only", status = "statistic_only"
  )
}

conditional_variable_test_rows <- function(conditional, data, calibration, B, adjust) {
  variables <- colnames(conditional)
  rows <- lapply(seq_along(variables), function(j) {
    rest <- setdiff(seq_along(variables), j)
    if (!length(rest)) return(conditional_table_template("variable"))
    result <- conditional_dependence_test(conditional[, j], data[, rest, drop = FALSE], calibration, B)
    data.frame(
      variable = variables[j], conditioning_variables = paste(variables[rest], collapse = ", "),
      representation = "conditional_standardized_residual_vs_original_conditioning_set",
      statistic = result$statistic, effect_size = result$effect_size, p_value = result$p_value,
      p_adjusted = NA_real_, method = result$method, calibration = result$calibration,
      null_hypothesis = paste0("R_", variables[j], " is independent of X_{-", variables[j], "}."),
      status = result$status, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (nrow(out)) out$p_adjusted <- stats::p.adjust(out$p_value, method = adjust)
  out
}

conditional_pair_test_rows <- function(conditional, data, calibration, B, adjust) {
  variables <- colnames(conditional)
  rows <- lapply(seq_along(variables), function(j) {
    lapply(setdiff(seq_along(variables), j), function(k) {
      result <- conditional_dependence_test(conditional[, j], data[, k, drop = FALSE], calibration, B)
      data.frame(
        residual_variable = variables[j], predictor = variables[k],
        representation = "conditional_standardized_residual_vs_original_predictor",
        statistic = result$statistic, effect_size = result$effect_size, p_value = result$p_value,
        p_adjusted = NA_real_, method = result$method, calibration = result$calibration,
        null_hypothesis = paste0("R_", variables[j], " is independent of X_", variables[k], "."),
        status = result$status, stringsAsFactors = FALSE
      )
    })
  })
  out <- do.call(rbind, unlist(rows, recursive = FALSE))
  if (nrow(out)) out$p_adjusted <- stats::p.adjust(out$p_value, method = adjust)
  out
}

conditional_scale_test_rows <- function(conditional, data, calibration, B, adjust) {
  variables <- colnames(conditional)
  rows <- lapply(seq_along(variables), function(j) {
    rest <- setdiff(seq_along(variables), j)
    if (!length(rest)) return(conditional_table_template("scale"))
    result <- conditional_dependence_test(conditional[, j]^2, data[, rest, drop = FALSE], calibration, B)
    data.frame(
      variable = variables[j], conditioning_variables = paste(variables[rest], collapse = ", "),
      representation = "squared_conditional_standardized_residual_vs_original_conditioning_set",
      statistic = result$statistic, effect_size = result$effect_size, p_value = result$p_value,
      p_adjusted = NA_real_, method = result$method, calibration = result$calibration,
      null_hypothesis = paste0("R_", variables[j], "^2 is independent of X_{-", variables[j], "}."),
      status = result$status, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (nrow(out)) out$p_adjusted <- stats::p.adjust(out$p_value, method = adjust)
  out
}

pairwise_independence_rows <- function(whitened, calibration, B, adjust) {
  variables <- colnames(whitened)
  pairs <- utils::combn(seq_along(variables), 2L)
  rows <- lapply(seq_len(ncol(pairs)), function(k) {
    i <- pairs[1L, k]
    j <- pairs[2L, k]
    if (requireNamespace("energy", quietly = TRUE)) {
      if (calibration %in% c("parametric_fixed", "parametric_refit")) {
        statistic <- unname(energy::dcov(whitened[, i], whitened[, j]))
        p_value <- NA_real_
        method <- "energy::dcov (parametric-bootstrap rank)"
        status <- "bootstrap_rank_pending"
      } else {
        test <- energy::dcov.test(whitened[, i], whitened[, j], R = B)
        statistic <- unname(test$statistic)
        p_value <- test$p.value
        method <- "energy::dcov.test"
        status <- "ok"
      }
    } else {
      statistic <- stats::cor(whitened[, i], whitened[, j], method = "spearman")
      p_value <- NA_real_
      method <- "Spearman correlation proxy; install energy for distance covariance"
      status <- "statistic_only"
    }
    data.frame(
      var1 = variables[i],
      var2 = variables[j],
      representation = "symmetric_zca_coordinate_pair",
      statistic = statistic,
      p_value = p_value,
      method = method,
      status = status
    )
  })
  out <- do.call(rbind, rows)
  out$p_adjusted <- stats::p.adjust(out$p_value, method = adjust)
  out
}

pairwise_global_row <- function(pair_rows) {
  if (is.null(pair_rows) || !nrow(pair_rows)) {
    return(diagnostic_row(
      test = "pairwise_independence",
      target = "full_gaussian_observed_data_fit",
      statistic = NA_real_,
      p_value = NA_real_,
      method = "no pairs available",
      status = "unavailable",
      representation = "symmetric_zca"
    ))
  }
  diagnostic_row(
    test = "pairwise_independence",
    target = "full_gaussian_observed_data_fit",
    statistic = max_or_na(pair_rows$statistic),
    p_value = min_or_na(pair_rows$p_adjusted),
    method = "minimum adjusted pairwise p-value",
    status = unique(pair_rows$status)[1],
    representation = "symmetric_zca"
  )
}

diagnostic_row <- function(test, target, statistic, p_value, method, status,
                           representation = "not_applicable") {
  data.frame(
    test = test,
    target = target,
    representation = representation,
    statistic = as.numeric(statistic),
    p_value = as.numeric(p_value),
    method = method,
    status = status,
    stringsAsFactors = FALSE
  )
}

sem_misfit_targets <- function() {
  c(
    mean_covariance_fit = "Model-implied observed means and covariance reproduce the first two moments.",
    full_gaussian_observed_data_fit = "Observed data follow the fitted multivariate Gaussian mean/covariance model.",
    conditional_adequacy_named_variable = "A named variable follows its fitted conditional distribution given named conditioning variables."
  )
}

find_openmx_fitfunction_expectations <- function(object) {
  algebras <- object$output$algebras
  for (nm in names(algebras)) {
    attrs <- attributes(algebras[[nm]])
    if (!is.null(attrs$expMean) && !is.null(attrs$expCov)) {
      return(list(expMean = attrs$expMean, expCov = attrs$expCov))
    }
  }
  stop("OpenMx expected mean/covariance attributes were not found.", call. = FALSE)
}

openmx_observed_data <- function(object) {
  if (!is.null(object$data$observed)) {
    return(object$data$observed)
  }
  NULL
}

openmx_fit_summary <- function(object) {
  out <- object$output
  if (is.null(out)) return(NULL)
  c(
    chisq = out$Chi %||% NA_real_,
    df = out$degreesOfFreedom %||% NA_real_,
    pvalue = out$p %||% NA_real_
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

min_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

max_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}
