options(repos = c(CRAN = "https://cloud.r-project.org"))

required <- c("lavaan", "psychTools", "energy", "dHSIC", "goftest", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

source("sem_misfit_api_prototype.R")

screen_B <- 19
set.seed(7707)

make_complete_numeric <- function(data, vars) {
  out <- as.data.frame(data[, vars, drop = FALSE])
  out <- out[stats::complete.cases(out), , drop = FALSE]
  out[] <- lapply(out, as.numeric)
  out
}

fit_screen <- function(name, source, data, model, notes) {
  fit <- lavaan::cfa(
    model,
    data = data,
    meanstructure = TRUE,
    std.lv = TRUE,
    estimator = "ML",
    missing = "listwise"
  )
  diag <- sem_misfit(
    fit,
    tests = "default",
    calibration = "permutation",
    B = screen_B,
    adjust = "holm",
    seed = 7707
  )
  tests <- as.data.frame(diag)
  fitmeas <- lavaan::fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"))
  top_pairs <- head(diag$pairs[order(diag$pairs$p_adjusted, -diag$pairs$statistic), ], 5)
  top_vars <- head(diag$variables[order(diag$variables$p_adjusted, -diag$variables$statistic), ], 5)

  list(
    summary = data.frame(
      candidate = name,
      source = source,
      n_available = attr(data, "n_available") %||% nrow(data),
      n_analyzed = nrow(data),
      p = ncol(data),
      chisq = unname(fitmeas["chisq"]),
      df = unname(fitmeas["df"]),
      sem_p = unname(fitmeas["pvalue"]),
      cfi = unname(fitmeas["cfi"]),
      rmsea = unname(fitmeas["rmsea"]),
      srmr = unname(fitmeas["srmr"]),
      global_distribution_p = tests$p_value[tests$test == "global_distribution"],
      global_independence_p = tests$p_value[tests$test == "global_independence"],
      pairwise_min_padj = tests$p_value[tests$test == "pairwise_independence"],
      marginal_min_padj = tests$p_value[tests$test == "marginal_distribution"],
      row_distance_p = tests$p_value[tests$test == "case_distance"],
      notes = notes,
      stringsAsFactors = FALSE
    ),
    tests = transform(tests, candidate = name),
    pairs = transform(top_pairs, candidate = name),
    variables = transform(top_vars, candidate = name),
    fit = fit,
    diag = diag
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

library(psychTools)

data(bfi, package = "psychTools")
bfi_vars <- c(
  "A2", "A3", "A4", "A5",
  "C1", "C2", "C3", "C4",
  "E1", "E2", "E3", "E4",
  "N1", "N2", "N3", "N4",
  "O1", "O2", "O3", "O4"
)
bfi_data <- make_complete_numeric(bfi, bfi_vars)
attr(bfi_data, "n_available") <- nrow(bfi_data)
bfi_model <- "
Agree =~ A2 + A3 + A4 + A5
Conscientious =~ C1 + C2 + C3 + C4
Extraversion =~ E1 + E2 + E3 + E4
Neuroticism =~ N1 + N2 + N3 + N4
Openness =~ O1 + O2 + O3 + O4
"

data(spi, package = "psychTools")
spi_vars <- c(
  "q_90", "q_1763", "q_253",
  "q_1290", "q_1744", "q_1979",
  "q_979", "q_4252", "q_1989",
  "q_1904", "q_4243", "q_1416",
  "q_128", "q_2745", "q_2754"
)
spi_data <- make_complete_numeric(spi, spi_vars)
attr(spi_data, "n_available") <- nrow(spi_data)
spi_model <- "
Agree =~ q_90 + q_1763 + q_253
Conscientious =~ q_1290 + q_1744 + q_1979
Neuroticism =~ q_979 + q_4252 + q_1989
Extraversion =~ q_1904 + q_4243 + q_1416
Openness =~ q_128 + q_2745 + q_2754
"

data(msq, package = "psychTools")
msq_vars <- c(
  "active", "energetic", "vigorous", "full.of.pep",
  "happy", "cheerful", "delighted", "pleased",
  "afraid", "nervous", "scared", "ashamed",
  "angry", "hostile", "grouchy", "frustrated"
)
msq_data <- make_complete_numeric(msq, msq_vars)
attr(msq_data, "n_available") <- nrow(msq_data)
msq_model <- "
Energy =~ active + energetic + vigorous + full.of.pep
PositiveAffect =~ happy + cheerful + delighted + pleased
FearGuilt =~ afraid + nervous + scared + ashamed
Hostility =~ angry + hostile + grouchy + frustrated
"

data(sai, package = "psychTools")
sai_vars <- c(
  "calm", "secure", "at.ease", "rested", "comfortable", "confident", "relaxed", "content",
  "tense", "upset", "worrying", "anxious", "nervous", "jittery", "high.strung", "worried"
)
sai_time1 <- sai[sai$time == 1, , drop = FALSE]
sai_data <- make_complete_numeric(sai_time1, sai_vars)
attr(sai_data, "n_available") <- nrow(sai_data)
sai_model <- "
CalmPositive =~ calm + secure + at.ease + rested + comfortable + confident + relaxed + content
AnxietyNegative =~ tense + upset + worrying + anxious + nervous + jittery + high.strung + worried
"

data(holzinger.swineford, package = "psychTools")
hs_vars <- c(
  "t01_visperc", "t02_cubes", "t03_frmbord",
  "t05_geninfo", "t06_paracomp", "t07_sentcomp",
  "t08_wordclas", "t09_wordmean",
  "t10_addition", "t12_countdot", "t13_sccaps"
)
hs_data <- make_complete_numeric(holzinger.swineford, hs_vars)
attr(hs_data, "n_available") <- nrow(hs_data)
hs_model <- "
Visual =~ t01_visperc + t02_cubes + t03_frmbord
Verbal =~ t05_geninfo + t06_paracomp + t07_sentcomp + t08_wordclas + t09_wordmean
Speed =~ t10_addition + t12_countdot + t13_sccaps
"

candidates <- list(
  BFI = list(
    source = "psychTools::bfi",
    data = bfi_data,
    model = bfi_model,
    notes = "Open Big Five item data; large, familiar, Likert-as-continuous caveat."
  ),
  SPI = list(
    source = "psychTools::spi",
    data = spi_data,
    model = spi_model,
    notes = "Open SAPA personality items; largest and modern, less familiar to SEM readers."
  ),
  MSQ = list(
    source = "psychTools::msq",
    data = msq_data,
    model = msq_model,
    notes = "Open affect-item data; strong substantive nonlinearity/subgroup possibilities."
  ),
  SAI = list(
    source = "psychTools::sai time == 1",
    data = sai_data,
    model = sai_model,
    notes = "Open single-occasion state-anxiety items; large and interpretable, but floor effects likely dominate."
  ),
  HolzingerSwineford = list(
    source = "psychTools::holzinger.swineford",
    data = hs_data,
    model = hs_model,
    notes = "Classic SEM dataset; familiar but small for nonparametric residual diagnostics."
  )
)

results <- lapply(names(candidates), function(name) {
  x <- candidates[[name]]
  message("Screening ", name)
  fit_screen(name, x$source, x$data, x$model, x$notes)
})
names(results) <- names(candidates)

summary_df <- do.call(rbind, lapply(results, `[[`, "summary"))
tests_df <- do.call(rbind, lapply(results, `[[`, "tests"))
pairs_df <- do.call(rbind, lapply(results, `[[`, "pairs"))
variables_df <- do.call(rbind, lapply(results, `[[`, "variables"))

write.csv(summary_df, "empirical_example_candidates.csv", row.names = FALSE)
write.csv(tests_df, "empirical_example_candidate_tests.csv", row.names = FALSE)
write.csv(pairs_df, "empirical_example_candidate_top_pairs.csv", row.names = FALSE)
write.csv(variables_df, "empirical_example_candidate_top_variables.csv", row.names = FALSE)

print(summary_df)
print(pairs_df)
print(variables_df)
