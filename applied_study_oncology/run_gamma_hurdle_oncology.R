# Gamma-hurdle benchmark for the oncology applied study.
#
# The analytic cohort, outcome, treatment coding, baseline adjustment set, and
# subgroup definitions match run_zicbcf_oncology.R.  The model is the
# gamma_hurdle benchmark implemented in zicbcf: a Bernoulli model for Y > 0
# and a Gamma model for the positive outcome.  It returns posterior draws of
# the response-scale treatment effect as well as p0 and p1 for the hurdle
# (participation) margin.

library(zicbcf)

outdir <- "applied_study_oncology"
source(file.path(outdir, "oncology_common.R"))

N_CHAINS <- 4L
N_BURN <- 1000L
N_SIM <- 1000L
CHAIN_SEEDS <- seq_len(N_CHAINS)

data_file <- file.path(outdir, "zic_bcf_headneck_analysis_data.csv")
if (!file.exists(data_file)) stop("Missing analysis data: ", data_file)

cat("Loading dataset...\n")
df <- read.csv(data_file, stringsAsFactors = FALSE)
y <- df$cumulative_severe_ae_duration
z <- as.integer(df$treatment)

if (anyNA(y) || anyNA(z) || any(y < 0) || !all(z %in% c(0L, 1L))) {
  stop("Outcome must be nonnegative and treatment must be coded 0/1 without missing values.")
}

# Dummy encoding exactly mirrors the ZIC-BCF-Smear analysis, drawn from the
# same shared definition so that the two specifications cannot drift apart.
X <- onc_design_matrix(df)

cat(sprintf("Sample size: %d; zero outcomes: %d (%.1f%%)\n",
            nrow(df), sum(y == 0), 100 * mean(y == 0)))
cat("Design matrix columns:", paste(colnames(X), collapse = ", "), "\n")
cat(sprintf("Fitting Gamma Hurdle benchmark: %d chains of %d retained draws...\n",
            N_CHAINS, N_SIM))

# Settings match the Gamma-hurdle run in applied_study/run_gamma_hurdle.R, now
# replicated across chains so that convergence can be assessed rather than
# assumed.
chains <- lapply(CHAIN_SEEDS, function(seed) {
  set.seed(seed)
  cat("  chain", seed, "\n")
  gamma_hurdle(y = y, z = z, x = X, nburn = N_BURN, nsim = N_SIM, nthin = 1)
})

convergence <- zicbcf_convergence_table(list(
  `ATE (response scale)` = lapply(chains, function(f) f$ate),
  `Hurdle ATE (probability)` = lapply(chains, function(f) rowMeans(f$p1 - f$p0)),
  `Mean CATE across units` = lapply(chains, function(f) rowMeans(f$cate))
))
write.csv(convergence, file.path(outdir, "gamma_hurdle_convergence_diagnostics.csv"),
          row.names = FALSE)
cat(sprintf("  max Rhat = %.4f, min ESS = %.0f, max |Geweke z| = %.2f\n",
            max(convergence$rhat, na.rm = TRUE), min(convergence$ess, na.rm = TRUE),
            max(convergence$geweke_z, na.rm = TRUE)))

cate_draws <- do.call(rbind, lapply(chains, function(f) f$cate))
ate_draws <- unlist(lapply(chains, function(f) f$ate), use.names = FALSE)
hurdle_cate_draws <- do.call(rbind, lapply(chains, function(f) f$p1 - f$p0))
hurdle_ate_draws <- rowMeans(hurdle_cate_draws)

if (ncol(cate_draws) != nrow(df) || ncol(hurdle_cate_draws) != nrow(df)) {
  stop("Unexpected posterior-draw dimensions returned by gamma_hurdle().")
}

cate_mean <- colMeans(cate_draws)
hurdle_cate_mean <- colMeans(hurdle_cate_draws)
ate_ci <- quantile(ate_draws, c(0.025, 0.975))
hurdle_ate_ci <- quantile(hurdle_ate_draws, c(0.025, 0.975))

res <- data.frame(
  Metric = c(
    "ATE", "ATE Lower 95% CI", "ATE Upper 95% CI", "P(ATE>0)",
    "Mean CATE", "Median CATE",
    "Hurdle ATE", "Hurdle ATE Lower 95% CI", "Hurdle ATE Upper 95% CI",
    "P(Hurdle ATE>0)", "Mean Hurdle CATE", "Median Hurdle CATE"
  ),
  Value = c(
    mean(ate_draws), ate_ci[1], ate_ci[2], mean(ate_draws > 0),
    mean(cate_mean), median(cate_mean),
    mean(hurdle_ate_draws), hurdle_ate_ci[1], hurdle_ate_ci[2],
    mean(hurdle_ate_draws > 0), mean(hurdle_cate_mean), median(hurdle_cate_mean)
  )
)
write.csv(res, file.path(outdir, "gamma_hurdle_ate_results.csv"), row.names = FALSE)

saveRDS(
  list(
    cate_draws = cate_draws,
    hurdle_cate_draws = hurdle_cate_draws,
    ate_draws = ate_draws,
    hurdle_ate_draws = hurdle_ate_draws,
    cate_mean = cate_mean,
    hurdle_cate_mean = hurdle_cate_mean,
    model_settings = list(nburn = N_BURN, nsim = N_SIM, nthin = 1,
                          n_chains = N_CHAINS),
    seed = CHAIN_SEEDS
  ),
  file.path(outdir, "gamma_hurdle_draws.rds")
)

png(file.path(outdir, "18_ONC_gamma_hurdle_cate_histogram.png"), width = 900, height = 600)
hist(cate_mean, breaks = 50, col = "skyblue", border = "white",
     main = "Gamma Hurdle: unit-level CATE distribution",
     xlab = "CATE: effect of panitumumab on cumulative grade 3+ AE duration (days)")
abline(v = mean(ate_draws), col = "firebrick", lwd = 2, lty = 3)
dev.off()

png(file.path(outdir, "19_ONC_gamma_hurdle_hurdle_cate_histogram.png"), width = 900, height = 600)
hist(hurdle_cate_mean, breaks = 50, col = "lightgreen", border = "white",
     main = "Gamma Hurdle: unit-level hurdle-margin CATE distribution",
     xlab = "CATE: effect of panitumumab on Pr(any grade 3+ AE)")
abline(v = mean(hurdle_ate_draws), col = "firebrick", lwd = 2, lty = 3)
dev.off()

cat(sprintf("ATE: %.3f days, 95%% CrI [%.3f, %.3f], P(ATE>0)=%.3f\n",
            mean(ate_draws), ate_ci[1], ate_ci[2], mean(ate_draws > 0)))
cat(sprintf("Hurdle ATE: %.4f, 95%% CrI [%.4f, %.4f], P(>0)=%.3f\n",
            mean(hurdle_ate_draws), hurdle_ate_ci[1], hurdle_ate_ci[2],
            mean(hurdle_ate_draws > 0)))
cat("Completed successfully.\n")
