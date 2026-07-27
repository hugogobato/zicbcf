# Matched ZIC-BCF-Smear and Gamma-hurdle analysis of distinct on-study days
# with at least one grade 3+ adverse event.  Run build_nonoverlap_calendar_day_endpoint.R first.

suppressPackageStartupMessages({
  library(zicbcf)
  library(dplyr)
  library(readr)
})

script_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))])
if (length(script_file) != 1L) stop("Run this script with Rscript.")

data_dir <- dirname(normalizePath(script_file))
out_dir <- file.path(data_dir, "analysis")
data_file <- file.path(out_dir, "nonoverlap_calendar_day_analysis_data.csv")
if (!file.exists(data_file)) stop("Run build_nonoverlap_calendar_day_endpoint.R first.")

df <- read_csv(data_file, show_col_types = FALSE)
y <- df$nonoverlap_severe_ae_days
z <- as.integer(df$treatment)

if (anyNA(y) || any(y < 0) || anyNA(z) || !all(z %in% c(0L, 1L))) {
  stop("Outcome must be nonnegative and treatment must be complete 0/1 coding.")
}

X <- model.matrix(
  ~ age + b_ecogct + sex + diagtype + tumcat + hpv + dstatus,
  data = df
)[, -1, drop = FALSE]
storage.mode(X) <- "double"

if (anyNA(X)) stop("Design matrix contains missing values.")

subgroup_labels <- list(
  Sex = factor(df$sex, levels = c("Male", "Female")),
  `ECOG status` = factor(df$b_ecogct, levels = c(0, 1), labels = c("ECOG 0", "ECOG 1")),
  `Diagnosis site` = factor(df$diagtype),
  `Tumor category` = factor(df$tumcat, levels = c("N", "Y"), labels = c("Tumor cat. N", "Tumor cat. Y")),
  `HPV status` = factor(df$hpv, levels = c("Negative", "Positive", "Unknown")),
  `Disease status` = factor(df$dstatus, levels = c("newly diagnosed", "recurrent")),
  `Age group` = cut(df$age, breaks = c(-Inf, 49, 59, 69, Inf), labels = c("<50", "50-59", "60-69", "70+"))
)

posterior_summary <- function(draws) {
  ci <- quantile(draws, c(0.025, 0.975))
  c(
    Estimate = mean(draws), CI_low = unname(ci[1L]), CI_high = unname(ci[2L]),
    P_gt_0 = mean(draws > 0)
  )
}

subgroup_table <- function(draws, model, estimand) {
  rows <- list(
    data.frame(
      Model = model, Estimand = estimand, Covariate = "Overall", Level = "All",
      N = ncol(draws), t(posterior_summary(rowMeans(draws))), check.names = FALSE
    )
  )

  for (group_name in names(subgroup_labels)) {
    label <- subgroup_labels[[group_name]]
    for (level in levels(label)) {
      index <- which(label == level)
      if (length(index) < 30L) next
      rows[[length(rows) + 1L]] <- data.frame(
        Model = model, Estimand = estimand, Covariate = group_name, Level = level,
        N = length(index), t(posterior_summary(rowMeans(draws[, index, drop = FALSE]))),
        check.names = FALSE
      )
    }
  }
  bind_rows(rows)
}

contrast_table <- function(draws, model, estimand) {
  rows <- list()
  for (group_name in names(subgroup_labels)) {
    label <- subgroup_labels[[group_name]]
    eligible_levels <- levels(label)[table(label) >= 30L]
    if (length(eligible_levels) < 2L) next
    comparisons <- combn(eligible_levels, 2L)
    for (column in seq_len(ncol(comparisons))) {
      first <- comparisons[1L, column]
      second <- comparisons[2L, column]
      contrast_draws <- rowMeans(draws[, which(label == first), drop = FALSE]) -
        rowMeans(draws[, which(label == second), drop = FALSE])
      rows[[length(rows) + 1L]] <- data.frame(
        Model = model, Estimand = estimand, Covariate = group_name,
        Contrast = paste(first, "-", second),
        t(posterior_summary(contrast_draws)), check.names = FALSE
      )
    }
  }
  bind_rows(rows)
}

cat(sprintf("Fitting models on %d patients, with %d zeros (%.1f%%).\n",
            nrow(df), sum(y == 0), 100 * mean(y == 0)))

# ZIC-BCF uses the same covariate set and posterior settings as the original
# record-duration analysis.  The estimated propensity is retained for the
# matched implementation, despite randomization.
set.seed(1)
ps_model <- glm(z ~ ., data = as.data.frame(X), family = binomial())
pihat <- predict(ps_model, type = "response")
zic_fit <- zicbcf_smear(
  y = y, z = z, x_control = X, x_moderate = X, pihat = pihat,
  nburn = 2000, nsim = 4000
)
zic_ate_draws <- zic_fit$ate
zic_cate_draws <- zic_fit$cate
zic_hurdle_cate_draws <- pnorm(zic_fit$mu_b + zic_fit$tau_b) - pnorm(zic_fit$mu_b)
zic_hurdle_ate_draws <- rowMeans(zic_hurdle_cate_draws)

# The parametric Gamma-hurdle benchmark matches the original benchmark's
# covariates, burn-in, posterior draws, thinning, and random seed.
set.seed(1)
gamma_fit <- gamma_hurdle(y = y, z = z, x = X, nburn = 1000, nsim = 1000, nthin = 1)
gamma_ate_draws <- gamma_fit$ate
gamma_cate_draws <- gamma_fit$cate
gamma_hurdle_cate_draws <- gamma_fit$p1 - gamma_fit$p0
gamma_hurdle_ate_draws <- rowMeans(gamma_hurdle_cate_draws)

overall_results <- bind_rows(
  data.frame(Model = "ZIC-BCF-Smear", Estimand = "Total non-overlapping severe-AE days", t(posterior_summary(zic_ate_draws))),
  data.frame(Model = "ZIC-BCF-Smear", Estimand = "Pr(any non-overlapping severe-AE day)", t(posterior_summary(zic_hurdle_ate_draws))),
  data.frame(Model = "Gamma hurdle", Estimand = "Total non-overlapping severe-AE days", t(posterior_summary(gamma_ate_draws))),
  data.frame(Model = "Gamma hurdle", Estimand = "Pr(any non-overlapping severe-AE day)", t(posterior_summary(gamma_hurdle_ate_draws)))
)

subgroups <- bind_rows(
  subgroup_table(zic_cate_draws, "ZIC-BCF-Smear", "Total non-overlapping severe-AE days"),
  subgroup_table(zic_hurdle_cate_draws, "ZIC-BCF-Smear", "Pr(any non-overlapping severe-AE day)"),
  subgroup_table(gamma_cate_draws, "Gamma hurdle", "Total non-overlapping severe-AE days"),
  subgroup_table(gamma_hurdle_cate_draws, "Gamma hurdle", "Pr(any non-overlapping severe-AE day)")
)

contrasts <- bind_rows(
  contrast_table(zic_cate_draws, "ZIC-BCF-Smear", "Total non-overlapping severe-AE days"),
  contrast_table(zic_hurdle_cate_draws, "ZIC-BCF-Smear", "Pr(any non-overlapping severe-AE day)"),
  contrast_table(gamma_cate_draws, "Gamma hurdle", "Total non-overlapping severe-AE days"),
  contrast_table(gamma_hurdle_cate_draws, "Gamma hurdle", "Pr(any non-overlapping severe-AE day)")
)

write_csv(overall_results, file.path(out_dir, "nonoverlap_overall_model_results.csv"))
write_csv(subgroups, file.path(out_dir, "nonoverlap_subgroup_model_results.csv"))
write_csv(contrasts, file.path(out_dir, "nonoverlap_subgroup_contrasts.csv"))

saveRDS(
  list(
    outcome = "Distinct on-study days with at least one grade 3+ AE",
    zic = list(ate = zic_ate_draws, cate = zic_cate_draws, hurdle_ate = zic_hurdle_ate_draws, hurdle_cate = zic_hurdle_cate_draws),
    gamma = list(ate = gamma_ate_draws, cate = gamma_cate_draws, hurdle_ate = gamma_hurdle_ate_draws, hurdle_cate = gamma_hurdle_cate_draws),
    zic_settings = list(nburn = 2000, nsim = 4000, seed = 1),
    gamma_settings = list(nburn = 1000, nsim = 1000, nthin = 1, seed = 1)
  ),
  file.path(out_dir, "nonoverlap_model_posterior_draws.rds")
)

png(file.path(out_dir, "nonoverlap_ate_posteriors.png"), width = 1000, height = 650)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
hist(zic_ate_draws, breaks = 40, col = "skyblue", border = "white",
     main = "ZIC-BCF-Smear", xlab = "ATE: distinct severe-AE days")
abline(v = c(0, mean(zic_ate_draws), quantile(zic_ate_draws, c(0.025, 0.975))),
       col = c("grey40", "firebrick", "firebrick"), lty = c(3, 1, 2), lwd = c(1, 2, 2))
hist(gamma_ate_draws, breaks = 40, col = "lightgreen", border = "white",
     main = "Gamma hurdle", xlab = "ATE: distinct severe-AE days")
abline(v = c(0, mean(gamma_ate_draws), quantile(gamma_ate_draws, c(0.025, 0.975))),
       col = c("grey40", "firebrick", "firebrick"), lty = c(3, 1, 2), lwd = c(1, 2, 2))
dev.off()

cat("Analysis completed.\n")
print(overall_results)
