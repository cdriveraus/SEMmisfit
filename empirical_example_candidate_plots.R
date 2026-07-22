required <- c("ISLR2", "NHANES", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

dir.create("empirical_example_figures", showWarnings = FALSE)

data(Wage, package = "ISLR2")
wage <- as.data.frame(Wage)
wage$education_num <- as.numeric(wage$education)
wage$linear_residual <- residuals(stats::lm(logwage ~ age + year + education_num, data = wage))

p_wage <- ggplot2::ggplot(wage, ggplot2::aes(age, linear_residual)) +
  ggplot2::geom_point(alpha = 0.18, size = 0.8) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 0.8) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  ggplot2::labs(
    x = "Age",
    y = "Linear-model residual for log wage",
    title = "Wage: residual curvature after linear path model"
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(
  "empirical_example_figures/wage_logwage_residual_by_age.png",
  p_wage,
  width = 6.5,
  height = 4.2,
  dpi = 180
)

data(NHANES, package = "NHANES")
nh <- as.data.frame(NHANES)
nh <- nh[nh$Age >= 18, ]
nh_vars <- c("Age", "BMI", "BPSysAve", "BPDiaAve", "DirectChol", "TotChol", "Pulse")
nh <- nh[, nh_vars]
nh <- nh[stats::complete.cases(nh), ]
nh$linear_residual <- residuals(stats::lm(
  BPSysAve ~ Age + BMI + BPDiaAve + DirectChol + TotChol + Pulse,
  data = nh
))

p_nhanes_age <- ggplot2::ggplot(nh, ggplot2::aes(Age, linear_residual)) +
  ggplot2::geom_point(alpha = 0.12, size = 0.65) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 0.8) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  ggplot2::labs(
    x = "Age",
    y = "Linear-model residual for systolic blood pressure",
    title = "NHANES: residual curvature by age"
  ) +
  ggplot2::theme_minimal(base_size = 11)

p_nhanes_bmi <- ggplot2::ggplot(nh, ggplot2::aes(BMI, linear_residual)) +
  ggplot2::geom_point(alpha = 0.12, size = 0.65) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 0.8) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  ggplot2::labs(
    x = "BMI",
    y = "Linear-model residual for systolic blood pressure",
    title = "NHANES: residual pattern by BMI"
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(
  "empirical_example_figures/nhanes_bpsys_residual_by_age.png",
  p_nhanes_age,
  width = 6.5,
  height = 4.2,
  dpi = 180
)

ggplot2::ggsave(
  "empirical_example_figures/nhanes_bpsys_residual_by_bmi.png",
  p_nhanes_bmi,
  width = 6.5,
  height = 4.2,
  dpi = 180
)

print(c(
  "empirical_example_figures/wage_logwage_residual_by_age.png",
  "empirical_example_figures/nhanes_bpsys_residual_by_age.png",
  "empirical_example_figures/nhanes_bpsys_residual_by_bmi.png"
))
