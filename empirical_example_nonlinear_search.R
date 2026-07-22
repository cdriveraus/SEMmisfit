options(repos = c(CRAN = "https://cloud.r-project.org"))

required <- c("lavaan", "energy", "dHSIC", "goftest", "mgcv", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

source("sem_misfit_api_prototype.R")

set.seed(770713)
pair_B <- 199
global_B <- 199
top_pair_count <- 5

make_fit <- function(model, data) {
  lavaan::sem(
    model,
    data = data,
    meanstructure = TRUE,
    estimator = "ML",
    missing = "listwise",
    fixed.x = FALSE
  )
}

pair_rank <- function(whitened) {
  pairs <- utils::combn(colnames(whitened), 2, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    x <- whitened[, pair[1]]
    y <- whitened[, pair[2]]
    data.frame(
      var1 = pair[1],
      var2 = pair[2],
      dcor = energy::dcor(x, y),
      pearson = stats::cor(x, y),
      spearman = stats::cor(x, y, method = "spearman"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(-out$dcor), , drop = FALSE]
}

test_top_pairs <- function(pair_stats, whitened, B = pair_B, top = top_pair_count) {
  top_rows <- head(pair_stats, top)
  tested <- lapply(seq_len(nrow(top_rows)), function(i) {
    row <- top_rows[i, ]
    test <- energy::dcov.test(whitened[, row$var1], whitened[, row$var2], R = B)
    data.frame(
      var1 = row$var1,
      var2 = row$var2,
      dcor = row$dcor,
      pearson = row$pearson,
      spearman = row$spearman,
      dcov_statistic = unname(test$statistic),
      dcov_p = test$p.value,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, tested)
  out$p_adjusted_top = stats::p.adjust(out$dcov_p, method = "holm")
  out
}

global_distribution <- function(whitened) {
  test <- energy::mvnorm.etest(whitened, R = global_B)
  c(statistic = unname(test$statistic), p = test$p.value)
}

global_independence <- function(whitened) {
  test <- dHSIC::dhsic.test(
    whitened,
    matrix.input = TRUE,
    B = global_B,
    method = "permutation"
  )
  c(statistic = unname(test$statistic), p = test$p.value)
}

row_distance_ad <- function(row_distance, df) {
  test <- goftest::ad.test(row_distance, "pchisq", df = df)
  c(statistic = unname(test$statistic), p = test$p.value)
}

marginal_ad <- function(whitened) {
  rows <- lapply(colnames(whitened), function(v) {
    test <- goftest::ad.test(whitened[, v], "pnorm", 0, 1)
    data.frame(variable = v, statistic = unname(test$statistic), p = test$p.value)
  })
  out <- do.call(rbind, rows)
  out$p_adjusted <- stats::p.adjust(out$p, method = "holm")
  out[order(out$p_adjusted, -out$statistic), , drop = FALSE]
}

smooth_check <- function(data, outcome, predictor, covariates = character()) {
  rhs <- paste(c(paste0("s(", predictor, ", k = 6)"), covariates), collapse = " + ")
  fit <- mgcv::gam(stats::as.formula(paste(outcome, "~", rhs)), data = data, method = "REML")
  st <- summary(fit)$s.table
  data.frame(
    outcome = outcome,
    predictor = predictor,
    edf = unname(st[1, "edf"]),
    p = unname(st[1, "p-value"]),
    stringsAsFactors = FALSE
  )
}

analyze_candidate <- function(name, source, data, model, note, smooth = NULL) {
  fit <- make_fit(model, data)
  expectation <- as_sem_expectation(fit)
  residuals <- compute_sem_residuals(expectation)
  whitened <- residuals$whitened

  pair_stats <- pair_rank(whitened)
  top_pairs <- test_top_pairs(pair_stats, whitened)
  gd <- global_distribution(whitened)
  gi <- global_independence(whitened)
  rd <- row_distance_ad(residuals$row_distance, ncol(whitened))
  mad <- marginal_ad(whitened)
  fitmeas <- lavaan::fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"))

  smooth_rows <- data.frame()
  if (!is.null(smooth)) {
    smooth_rows <- do.call(rbind, lapply(smooth, function(x) {
      out <- smooth_check(data, x$outcome, x$predictor, x$covariates %||% character())
      transform(out, candidate = name)
    }))
  }

  summary <- data.frame(
    candidate = name,
    source = source,
    n = nrow(data),
    p = ncol(data),
    chisq = unname(fitmeas["chisq"]),
    df = unname(fitmeas["df"]),
    sem_p = unname(fitmeas["pvalue"]),
    cfi = unname(fitmeas["cfi"]),
    rmsea = unname(fitmeas["rmsea"]),
    srmr = unname(fitmeas["srmr"]),
    global_distribution_p = gd["p"],
    global_independence_p = gi["p"],
    row_distance_p = rd["p"],
    marginal_min_padj = min(mad$p_adjusted, na.rm = TRUE),
    top_pair = paste(top_pairs$var1[1], top_pairs$var2[1], sep = "--"),
    top_pair_dcor = top_pairs$dcor[1],
    top_pair_p = top_pairs$dcov_p[1],
    note = note,
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    pairs = transform(top_pairs, candidate = name),
    marginal = transform(head(mad, 5), candidate = name),
    smooth = smooth_rows,
    fit = fit,
    residuals = residuals
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

candidates <- list()

if (requireNamespace("ISLR2", quietly = TRUE)) {
  data(Wage, package = "ISLR2")
  wage <- as.data.frame(Wage)
  wage$education_num <- as.numeric(wage$education)
  wage_data <- wage[, c("age", "year", "education_num", "logwage")]
  wage_data <- wage_data[stats::complete.cases(wage_data), ]
  candidates$WageLogPath <- list(
    source = "ISLR2::Wage",
    data = wage_data,
    model = "
      logwage ~ age + year + education_num
      age ~~ year + education_num
      year ~~ education_num
    ",
    note = "Observed-variable path model; saturated for linear mean/covariance, targeted at nonlinear age-wage residual dependence.",
    smooth = list(list(outcome = "logwage", predictor = "age", covariates = c("year", "education_num")))
  )
}

if (requireNamespace("NHANES", quietly = TRUE)) {
  data(NHANES, package = "NHANES")
  nh <- as.data.frame(NHANES)
  nh <- nh[nh$Age >= 18, ]
  nh_vars <- c("Age", "BMI", "BPSysAve", "BPDiaAve", "DirectChol", "TotChol", "Pulse")
  nh_data <- nh[, nh_vars]
  nh_data <- nh_data[stats::complete.cases(nh_data), ]
  candidates$NHANESBpPath <- list(
    source = "NHANES::NHANES adults",
    data = nh_data,
    model = "
      BPSysAve ~ Age + BMI + BPDiaAve + DirectChol + TotChol + Pulse
      Age ~~ BMI + BPDiaAve + DirectChol + TotChol + Pulse
      BMI ~~ BPDiaAve + DirectChol + TotChol + Pulse
      BPDiaAve ~~ DirectChol + TotChol + Pulse
      DirectChol ~~ TotChol + Pulse
      TotChol ~~ Pulse
    ",
    note = "Observed-variable health path model; saturated linear covariance baseline, targeted at nonlinear age/BMI/blood-pressure structure.",
    smooth = list(
      list(outcome = "BPSysAve", predictor = "Age", covariates = c("BMI", "BPDiaAve", "DirectChol", "TotChol", "Pulse")),
      list(outcome = "BPSysAve", predictor = "BMI", covariates = c("Age", "BPDiaAve", "DirectChol", "TotChol", "Pulse"))
    )
  )
}

if (requireNamespace("psychTools", quietly = TRUE)) {
  data(msq, package = "psychTools")
  msq_vars <- c(
    "active", "energetic", "vigorous", "full.of.pep",
    "happy", "cheerful", "delighted", "pleased",
    "afraid", "nervous", "scared", "ashamed",
    "angry", "hostile", "grouchy", "frustrated"
  )
  msq_data <- as.data.frame(msq[, msq_vars])
  msq_data <- msq_data[stats::complete.cases(msq_data), ]
  candidates$MSQAffectCfa <- list(
    source = "psychTools::msq",
    data = msq_data,
    model = "
      Energy =~ active + energetic + vigorous + full.of.pep
      PositiveAffect =~ happy + cheerful + delighted + pleased
      FearGuilt =~ afraid + nervous + scared + ashamed
      Hostility =~ angry + hostile + grouchy + frustrated
    ",
    note = "Affect-item CFA; likely local dependence and floor/skew structure rather than clean nonlinear cross-variable dependence."
  )

  data(spi, package = "psychTools")
  spi_vars <- c(
    "q_90", "q_1763", "q_253",
    "q_1290", "q_1744", "q_1979",
    "q_979", "q_4252", "q_1989",
    "q_1904", "q_4243", "q_1416",
    "q_128", "q_2745", "q_2754"
  )
  spi_data <- as.data.frame(spi[, spi_vars])
  spi_data <- spi_data[stats::complete.cases(spi_data), ]
  candidates$SPIBigFiveCfa <- list(
    source = "psychTools::spi",
    data = spi_data,
    model = "
      Agree =~ q_90 + q_1763 + q_253
      Conscientious =~ q_1290 + q_1744 + q_1979
      Neuroticism =~ q_979 + q_4252 + q_1989
      Extraversion =~ q_1904 + q_4243 + q_1416
      Openness =~ q_128 + q_2745 + q_2754
    ",
    note = "Personality-item CFA; strong N and fit, but top dependence may be same-facet local item dependence."
  )
}

results <- lapply(names(candidates), function(name) {
  message("Analyzing ", name)
  x <- candidates[[name]]
  analyze_candidate(name, x$source, x$data, x$model, x$note, x$smooth)
})
names(results) <- names(candidates)

summary_df <- do.call(rbind, lapply(results, `[[`, "summary"))
pairs_df <- do.call(rbind, lapply(results, `[[`, "pairs"))
marginal_df <- do.call(rbind, lapply(results, `[[`, "marginal"))
smooth_df <- do.call(rbind, lapply(results, `[[`, "smooth"))

write.csv(summary_df, "empirical_example_nonlinear_candidates.csv", row.names = FALSE)
write.csv(pairs_df, "empirical_example_nonlinear_top_pairs.csv", row.names = FALSE)
write.csv(marginal_df, "empirical_example_nonlinear_marginal.csv", row.names = FALSE)
write.csv(smooth_df, "empirical_example_nonlinear_smooth_checks.csv", row.names = FALSE)

print(summary_df)
print(pairs_df)
print(marginal_df)
print(smooth_df)
