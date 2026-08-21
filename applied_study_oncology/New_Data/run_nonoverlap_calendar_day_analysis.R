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
source(file.path(data_dir, "..", "oncology_common.R"))
out_dir <- file.path(data_dir, "analysis")
data_file <- file.path(out_dir, "nonoverlap_calendar_day_analysis_data.csv")
if (!file.exists(data_file)) stop("Run build_nonoverlap_calendar_day_endpoint.R first.")

df <- read_csv(data_file, show_col_types = FALSE)
y <- df$nonoverlap_severe_ae_days
z <- as.integer(df$treatment)

if (anyNA(y) || any(y < 0) || anyNA(z) || !all(z %in% c(0L, 1L))) {
  stop("Outcome must be nonnegative and treatment must be complete 0/1 coding.")
}

X <- onc_design_matrix(df)
subgroup_labels <- onc_subgroup_labels(df)

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
      if (length(index) < ONC_MIN_SUBGROUP_N) next
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
    eligible_levels <- levels(label)[table(label) >= ONC_MIN_SUBGROUP_N]
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

N_CHAINS <- 4L
CHAIN_SEEDS <- seq_len(N_CHAINS)

# ZIC-BCF uses the same covariate set and posterior settings as the
# record-duration analysis. The known randomized assignment probability is
# used directly rather than estimating a potentially separated propensity
# model from the realized arms.
pihat <- rep(0.5, length(z))
zic_chains <- lapply(CHAIN_SEEDS, function(seed) {
  set.seed(seed)
  cat("  ZIC-BCF-Smear chain", seed, "\n")
  zicbcf_smear(y = y, z = z, x_control = X, x_moderate = X, pihat = pihat,
               nburn = 2000, nsim = 4000)
})
zic_convergence <- zicbcf_convergence(zic_chains, n_cate_units = 10L)
zic_smearing <- zicbcf_smearing_diagnostics(zic_chains[[1L]], y = y, z = z, x = X)

zic_ate_draws <- unlist(lapply(zic_chains, function(f) f$ate), use.names = FALSE)
zic_cate_draws <- do.call(rbind, lapply(zic_chains, function(f) f$cate))
zic_hurdle_cate_draws <- do.call(rbind, lapply(zic_chains, function(f) {
  pnorm(f$mu_b + f$tau_b) - pnorm(f$mu_b)
}))
zic_hurdle_ate_draws <- rowMeans(zic_hurdle_cate_draws)

# The parametric Gamma-hurdle benchmark matches the original benchmark's
# covariates, burn-in, posterior draws, thinning, and chain configuration.
gamma_chains <- lapply(CHAIN_SEEDS, function(seed) {
  set.seed(seed)
  cat("  Gamma hurdle chain", seed, "\n")
  gamma_hurdle(y = y, z = z, x = X, nburn = 1000, nsim = 1000, nthin = 1)
})
gamma_convergence <- zicbcf_convergence_table(list(
  `ATE (response scale)` = lapply(gamma_chains, function(f) f$ate),
  `Hurdle ATE (probability)` = lapply(gamma_chains, function(f) rowMeans(f$p1 - f$p0)),
  `Mean CATE across units` = lapply(gamma_chains, function(f) rowMeans(f$cate))
))

gamma_ate_draws <- unlist(lapply(gamma_chains, function(f) f$ate), use.names = FALSE)
gamma_cate_draws <- do.call(rbind, lapply(gamma_chains, function(f) f$cate))
gamma_hurdle_cate_draws <- do.call(rbind, lapply(gamma_chains, function(f) f$p1 - f$p0))
gamma_hurdle_ate_draws <- rowMeans(gamma_hurdle_cate_draws)

write_csv(bind_rows(
  mutate(zic_convergence, Model = "ZIC-BCF-Smear"),
  mutate(gamma_convergence, Model = "Gamma hurdle")
), file.path(out_dir, "nonoverlap_convergence_diagnostics.csv"))
write_csv(bind_rows(
  mutate(zic_smearing$test, Quantity = "Breusch-Pagan on posterior-mean residuals"),
  mutate(zic_smearing$sensitivity, Quantity = "Locally smeared sensitivity")
), file.path(out_dir, "nonoverlap_smearing_diagnostics.csv"))
write_csv(zic_smearing$residual_scale,
          file.path(out_dir, "nonoverlap_smearing_residual_scale.csv"))

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
    zic_settings = list(nburn = 2000, nsim = 4000, seeds = CHAIN_SEEDS),
    gamma_settings = list(nburn = 1000, nsim = 1000, nthin = 1, seeds = CHAIN_SEEDS)
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
cat(sprintf("ZIC-BCF-Smear: max Rhat = %.4f, min ESS = %.0f, max |Geweke z| = %.2f\n",
            max(zic_convergence$rhat, na.rm = TRUE), min(zic_convergence$ess, na.rm = TRUE),
            max(zic_convergence$geweke_z, na.rm = TRUE)))
cat(sprintf("Gamma hurdle:  max Rhat = %.4f, min ESS = %.0f, max |Geweke z| = %.2f\n",
            max(gamma_convergence$rhat, na.rm = TRUE), min(gamma_convergence$ess, na.rm = TRUE),
            max(gamma_convergence$geweke_z, na.rm = TRUE)))
cat(sprintf("Smearing homoskedasticity: Breusch-Pagan p = %.4f; locally smeared ATE differs by %.1f%%\n",
            zic_smearing$test$p_value, 100 * zic_smearing$sensitivity$relative_change))
