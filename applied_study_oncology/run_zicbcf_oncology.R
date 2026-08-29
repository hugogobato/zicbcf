library(zicbcf)

# ==========================================================================
# Applied study (oncology): effect of RANDOMIZED ASSIGNMENT to panitumumab
# plus chemotherapy, versus chemotherapy alone, on the cumulative burden in
# days of severe (grade 3+) adverse events, in the SPECTRUM head-and-neck
# cancer trial.
#
# ESTIMAND. The target is the effect of assignment on accrued severe-adverse-
# event days under the trial's own differential treatment duration. Time on
# treatment is a post-randomization consequence of assignment (panitumumab
# continues as maintenance after chemotherapy stops), so it lies on the causal
# pathway from assignment to accrued burden. It is therefore deliberately
# absent from the adjustment set: conditioning on it would block part of the
# effect being estimated and would condition on a post-treatment variable. The
# consequence for interpretation is that the estimate is not a per-exposure
# rate and not a fixed-horizon restricted mean. The descriptive exposure
# analyses that quantify this are in
# New_Data/run_exposure_and_missingness_diagnostics.R and
# New_Data/run_exposure_rate_descriptive_analysis.R.
#
# The outcome is semicontinuous: a point mass at zero (patients with no
# grade 3+ AE, about 24%) plus a heavily right-skewed positive part in days,
# which is exactly the data structure ZIC-BCF-Smear targets.
#
# The trial is randomized, so the unadjusted arm difference is already
# unbiased for the average effect in this cohort. The propensity model and the
# baseline adjustment set are variance-reduction and precision devices, not
# requirements for validity.
#
# A subgroup CATE is obtained by averaging unit-level effects over the units in
# the subgroup *within each posterior draw*, which yields a posterior
# distribution for the subgroup effect and therefore a point estimate and a 95%
# credible interval. Between-group contrasts are computed the same way, so a
# contrast inherits the posterior covariance of its two subgroups.
# ==========================================================================

outdir <- "applied_study_oncology"
source(file.path(outdir, "oncology_common.R"))

N_CHAINS <- 4L
N_BURN <- 5000L
N_SIM <- 4000L
CHAIN_SEEDS <- seq_len(N_CHAINS)

cat("Loading dataset...\n")
data_file <- file.path(outdir, "zic_bcf_headneck_analysis_data.csv")
if (!file.exists(data_file)) stop("Run prepare_analysis_data.R first.")
df <- read.csv(data_file, stringsAsFactors = FALSE)
cat("Sample size:", nrow(df), "\n")

y <- df$cumulative_severe_ae_duration
z <- as.integer(df$treatment)

cat("Zeros in outcome:", sum(y == 0), "/", length(y),
    sprintf("(%.1f%%)\n", 100 * mean(y == 0)))

X <- onc_design_matrix(df)
cat("Design matrix columns:", paste(colnames(X), collapse = ", "), "\n")

# -- propensity score --------------------------------------------------------
# This is a randomized trial (260/260), so the design propensity is 0.5. Use
# that known assignment probability rather than estimating a potentially
# separated observational propensity model from the realized arms.
cat("Using the randomized design propensity pihat = 0.5...\n")
pihat <- rep(0.5, length(z))

# -- fit ZIC-BCF-Smear, one chain per seed -----------------------------------
cat(sprintf("Fitting ZIC-BCF-Smear: %d chains of %d retained draws...\n",
            N_CHAINS, N_SIM))
chains <- lapply(CHAIN_SEEDS, function(seed) {
  set.seed(seed)
  cat("  chain", seed, "\n")
  zicbcf_smear(y = y, z = z, x_control = X, x_moderate = X, pihat = pihat,
               nburn = N_BURN, nsim = N_SIM)
})

# -- convergence diagnostics -------------------------------------------------
cat("Computing convergence diagnostics...\n")
convergence <- zicbcf_convergence(chains, n_cate_units = 10L)
write.csv(convergence, file.path(outdir, "onc_convergence_diagnostics.csv"),
          row.names = FALSE)
cat(sprintf("  max Rhat = %.4f, min ESS = %.0f, max |Geweke z| = %.2f\n",
            max(convergence$rhat, na.rm = TRUE),
            min(convergence$ess, na.rm = TRUE),
            max(convergence$geweke_z, na.rm = TRUE)))

# -- smearing homoskedasticity diagnostics -----------------------------------
cat("Testing the homoskedasticity assumption behind Duan's smearing...\n")
smearing <- zicbcf_smearing_diagnostics(chains[[1L]], y = y, z = z, x = X)
write.csv(smearing$test, file.path(outdir, "onc_smearing_bp_test.csv"), row.names = FALSE)
write.csv(smearing$draw_test, file.path(outdir, "onc_smearing_bp_draws.csv"), row.names = FALSE)
write.csv(smearing$residual_scale, file.path(outdir, "onc_smearing_residual_scale.csv"),
          row.names = FALSE)
write.csv(smearing$sensitivity, file.path(outdir, "onc_smearing_sensitivity.csv"),
          row.names = FALSE)
cat(sprintf("  Breusch-Pagan p = %.4f; locally smeared ATE differs by %.1f%%\n",
            smearing$test$p_value, 100 * smearing$sensitivity$relative_change))

# -- pool the chains ---------------------------------------------------------
# The reported posterior is the pooled draw set across all chains.
cate_draws <- do.call(rbind, lapply(chains, function(f) f$cate))
ate_draws <- unlist(lapply(chains, function(f) f$ate), use.names = FALSE)
hurdle_cate_draws <- do.call(rbind, lapply(chains, function(f) {
  pnorm(f$mu_b + f$tau_b) - pnorm(f$mu_b)
}))
hurdle_ate_draws <- rowMeans(hurdle_cate_draws)

saveRDS(cate_draws, file.path(outdir, "onc_cate_draws.rds"))
saveRDS(list(hurdle_cate_draws = hurdle_cate_draws,
             hurdle_ate_draws = hurdle_ate_draws,
             hurdle_cate_mean = colMeans(hurdle_cate_draws)),
        file.path(outdir, "onc_hurdle_draws.rds"))

hurdle_ate_ci <- quantile(hurdle_ate_draws, c(0.025, 0.975))
write.csv(data.frame(
  Metric = c("Hurdle ATE", "Hurdle ATE Lower 95% CI",
             "Hurdle ATE Upper 95% CI", "P(Hurdle ATE>0)",
             "Mean Hurdle CATE", "Median Hurdle CATE"),
  Value = c(mean(hurdle_ate_draws), hurdle_ate_ci[1], hurdle_ate_ci[2],
            mean(hurdle_ate_draws > 0), mean(colMeans(hurdle_cate_draws)),
            median(colMeans(hurdle_cate_draws)))
), file.path(outdir, "onc_hurdle_ate_results.csv"), row.names = FALSE)

# -- ATE ---------------------------------------------------------------------
ate_mean <- mean(ate_draws)
ate_ci <- quantile(ate_draws, c(0.025, 0.975))
cat(sprintf("\nATE: %.2f days  95%% CrI [%.2f, %.2f]  P(>0)=%.3f\n",
            ate_mean, ate_ci[1], ate_ci[2], mean(ate_draws > 0)))

cate_mean <- colMeans(cate_draws)
write.csv(data.frame(
  Metric = c("ATE", "Lower 95% CI", "Upper 95% CI", "P(ATE>0)",
             "Mean CATE", "Median CATE"),
  Value = c(ate_mean, ate_ci[1], ate_ci[2], mean(ate_draws > 0),
            mean(cate_mean), median(cate_mean))
), file.path(outdir, "onc_ate_results.csv"), row.names = FALSE)

# -- unadjusted arm difference, for the randomization benchmark --------------
# Under randomization the unadjusted difference is itself unbiased for the
# cohort average effect, so proximity to it is the benchmark a model should
# meet, not a weakness of the model.
unadjusted <- mean(y[z == 1]) - mean(y[z == 0])
cat(sprintf("Unadjusted arm difference: %.2f days\n", unadjusted))
write.csv(data.frame(Metric = "Unadjusted arm difference (event-days)",
                     Value = unadjusted),
          file.path(outdir, "onc_unadjusted_difference.csv"), row.names = FALSE)

# ===========================================================================
# Subgroup CATEs and pairwise contrasts
# ===========================================================================
labels <- onc_subgroup_labels(df)

cate_table <- onc_subgroup_table(cate_draws, labels, "ZIC-BCF-Smear",
                                 "Summed severe-AE-record duration (event-days)")
names(cate_table)[names(cate_table) == "Estimate"] <- "CATE"
cate_table[, c("CATE", "CI_low", "CI_high")] <-
  round(cate_table[, c("CATE", "CI_low", "CI_high")], 2)
cate_table$P_gt_0 <- round(cate_table$P_gt_0, 3)
write.csv(cate_table, file.path(outdir, "onc_cate_subgroups.csv"), row.names = FALSE)
cat("\n=== Subgroup CATE estimates ===\n"); print(cate_table)

contrast_table <- onc_contrast_table(cate_draws, labels, "ZIC-BCF-Smear",
                                     "Summed severe-AE-record duration (event-days)")
names(contrast_table)[names(contrast_table) == "Estimate"] <- "Diff"
contrast_table[, c("Diff", "CI_low", "CI_high")] <-
  round(contrast_table[, c("Diff", "CI_low", "CI_high")], 2)
contrast_table$P_gt_0 <- round(contrast_table$P_gt_0, 3)
write.csv(contrast_table, file.path(outdir, "onc_cate_contrasts.csv"), row.names = FALSE)
cat("\n=== Between-group CATE contrasts (heterogeneity) ===\n"); print(contrast_table)

# ===========================================================================
# Figures
# ===========================================================================
XLAB <- "CATE: effect of panitumumab on cumulative grade 3+ AE duration (days)"

plot_df <- cate_table[cate_table$Covariate != "Overall", ]
plot_df$row_lab <- paste0(plot_df$Covariate, ": ", plot_df$Level)
plot_df <- plot_df[nrow(plot_df):1, ]
overall_cate <- cate_table$CATE[cate_table$Covariate == "Overall"]

png(file.path(outdir, "11_ONC_cate_subgroups.png"), width = 1000, height = 1100)
par(mar = c(5, 16, 4, 2))
yy <- seq_len(nrow(plot_df))
plot(plot_df$CATE, yy, xlim = range(c(plot_df$CI_low, plot_df$CI_high, 0)),
     yaxt = "n", pch = 19, col = "steelblue", cex = 1.3,
     xlab = XLAB, ylab = "",
     main = "Subgroup CATE estimates with 95% credible intervals")
axis(2, at = yy, labels = plot_df$row_lab, las = 1, cex.axis = 0.9)
segments(plot_df$CI_low, yy, plot_df$CI_high, yy, col = "steelblue", lwd = 2)
abline(v = 0, lty = 2, col = "grey50")
abline(v = overall_cate, lty = 3, col = "firebrick", lwd = 2)
legend("bottomright",
       legend = c("Subgroup CATE (95% CrI)", "No effect",
                  sprintf("Overall ATE = %.1f days", overall_cate)),
       col = c("steelblue", "grey50", "firebrick"),
       pch = c(19, NA, NA), lty = c(NA, 2, 3), lwd = c(2, 1, 2), bty = "n")
dev.off()

png(file.path(outdir, "12_ONC_ate_posterior.png"), width = 900, height = 600)
hist(ate_draws, breaks = 40, col = "skyblue", border = "white",
     main = "Posterior distribution of the ATE",
     xlab = "ATE: effect of panitumumab on cumulative grade 3+ AE duration (days)")
abline(v = ate_mean, col = "firebrick", lwd = 2)
abline(v = ate_ci, col = "firebrick", lwd = 2, lty = 2)
abline(v = 0, col = "grey40", lwd = 1, lty = 3)
legend("topright",
       legend = c(sprintf("Posterior mean = %.1f", ate_mean),
                  sprintf("95%% CrI [%.1f, %.1f]", ate_ci[1], ate_ci[2])),
       col = "firebrick", lwd = 2, lty = c(1, 2), bty = "n")
dev.off()

png(file.path(outdir, "13_ONC_cate_unit_histogram.png"), width = 900, height = 600)
hist(cate_mean, breaks = 50, col = "skyblue", border = "white",
     main = "Distribution of unit-level Conditional Average Treatment Effects",
     xlab = XLAB)
abline(v = ate_mean, col = "firebrick", lwd = 2, lty = 3)
dev.off()

dens_groups <- labels[c("ECOG status", "Disease status")]
png(file.path(outdir, "14_ONC_cate_subgroup_densities.png"), width = 1000, height = 850)
par(mfrow = c(2, 2), mar = c(4.5, 4, 3, 1))
for (gname in names(dens_groups)) {
  lab <- dens_groups[[gname]]
  lvs <- levels(lab)[table(lab) >= ONC_MIN_SUBGROUP_N]
  dens <- lapply(lvs, function(lv) density(rowMeans(cate_draws[, which(lab == lv), drop = FALSE])))
  xr <- range(sapply(dens, function(d) range(d$x)))
  yr <- range(sapply(dens, function(d) range(d$y)))
  cols <- c("steelblue", "firebrick", "darkgreen")[seq_along(lvs)]
  plot(NA, xlim = xr, ylim = yr, xlab = "Subgroup CATE (days)",
       ylab = "Posterior density", main = gname)
  for (i in seq_along(dens)) {
    lines(dens[[i]], col = cols[i], lwd = 2)
    polygon(dens[[i]], col = adjustcolor(cols[i], 0.2), border = NA)
  }
  abline(v = 0, lty = 3, col = "grey50")
  legend("topright", legend = lvs, col = cols, lwd = 2, bty = "n")
}
dev.off()

png(file.path(outdir, "15_ONC_cate_vs_age.png"), width = 900, height = 600)
plot(df$age, cate_mean, pch = 19, col = adjustcolor("steelblue", 0.4),
     xlab = "Baseline age (years)", ylab = "Unit posterior-mean CATE (days)",
     main = "Unit-level CATE versus baseline age")
lines(lowess(df$age, cate_mean, f = 0.5), col = "firebrick", lwd = 3)
abline(h = ate_mean, lty = 3, col = "grey40")
dev.off()

ct <- contrast_table
ct$row_lab <- paste0(ct$Covariate, ": ", ct$Contrast)
ct <- ct[nrow(ct):1, ]
credible <- ct$CI_low > 0 | ct$CI_high < 0
png(file.path(outdir, "16_ONC_cate_contrasts_forest.png"), width = 1000, height = 1100)
par(mar = c(5, 16, 4, 2))
yy <- seq_len(nrow(ct))
cols <- ifelse(credible, "firebrick", "grey50")
plot(ct$Diff, yy, xlim = range(c(ct$CI_low, ct$CI_high, 0)),
     yaxt = "n", pch = 19, col = cols, cex = 1.2,
     xlab = "Between-subgroup difference in CATE (days)", ylab = "",
     main = "Pairwise subgroup CATE contrasts with 95% credible intervals")
axis(2, at = yy, labels = ct$row_lab, las = 1, cex.axis = 0.75)
segments(ct$CI_low, yy, ct$CI_high, yy, col = cols, lwd = 2)
abline(v = 0, lty = 2, col = "grey40")
legend("bottomright",
       legend = c("Credible difference (CrI excludes 0)", "Not credible"),
       col = c("firebrick", "grey50"), pch = 19, bty = "n")
dev.off()

cat("\nATE:", round(ate_mean, 2),
    " 95% CrI [", round(ate_ci[1], 2), ",", round(ate_ci[2], 2), "]\n")
cat("Completed successfully.\n")
