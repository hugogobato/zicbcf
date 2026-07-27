library(zicbcf)

# ==========================================================================
# Applied study (oncology): effect of adding panitumumab to chemotherapy on
# the cumulative burden (in days) of severe (grade 3+) adverse events in the
# HeadNe_Amgen_2007_265 head-and-neck cancer trial.
#
# The outcome (cumulative_severe_ae_duration) is semicontinuous: a point mass
# at zero (patients with no grade 3+ AE, ~24%) plus a heavily right-skewed
# positive part (days), which is exactly the data structure ZIC-BCF targets.
#
# As in the MEPS dental study, the ZIC-BCF-Smear fit returns a full posterior
# of unit-level CATEs (cate_draws: nsim x n). A subgroup CATE is obtained by
# averaging the unit-level effects over the units in the subgroup *within each
# posterior draw*, which yields a posterior distribution for the subgroup
# effect and therefore a point estimate AND a 95% credible interval.
# Between-group contrasts are computed the same way, giving a direct posterior
# test for treatment-effect heterogeneity.
# ==========================================================================

set.seed(1)
outdir <- "applied_study_oncology"

cat("Loading dataset...\n")
df <- read.csv(file.path(outdir, "zic_bcf_headneck_analysis_data.csv"),
               stringsAsFactors = FALSE)
cat("Sample size:", nrow(df), "\n")

# -- outcome, treatment, covariates -----------------------------------------
y <- df$cumulative_severe_ae_duration           # cumulative grade 3+ AE days
z <- as.integer(df$treatment)                   # 1 = panitumumab + chemo, 0 = chemo

cat("Zeros in outcome:", sum(y == 0), "/", length(y),
    sprintf("(%.1f%%)\n", 100 * mean(y == 0)))

# Build a numeric design matrix (dummy-encode nominal factors; keep age and
# ECOG numeric). age and b_ecogct are numeric; the rest are categorical.
X <- model.matrix(
  ~ age + b_ecogct + sex + diagtype + tumcat + hpv + dstatus,
  data = df
)[, -1, drop = FALSE]                           # drop intercept column
storage.mode(X) <- "double"
cat("Design matrix columns:", paste(colnames(X), collapse = ", "), "\n")

# -- propensity score --------------------------------------------------------
# This is a randomized trial (260/260), so the propensity is ~0.5 by design;
# we still estimate it from covariates for consistency with the MEPS pipeline.
cat("Estimating propensity scores...\n")
ps_model <- glm(z ~ ., data = as.data.frame(X), family = binomial())
pihat <- predict(ps_model, type = "response")
cat("pihat range: [", round(min(pihat), 3), ",", round(max(pihat), 3), "]\n")

# -- fit ZIC-BCF-Smear -------------------------------------------------------
cat("Fitting ZIC-BCF-Smear...\n")
fit <- zicbcf_smear(
  y = y, z = z,
  x_control = X, x_moderate = X,
  pihat = pihat,
  nburn = 2000, nsim = 4000
)

cate_draws <- fit$cate            # nsim x n posterior of unit-level CATE
ate_draws  <- fit$ate
saveRDS(cate_draws, file.path(outdir, "onc_cate_draws.rds"))

# Save the participation (hurdle) posterior as well.  The fitted probit
# hurdle supplies p_0(x)=Phi(mu_b(x)) and p_1(x)=Phi(mu_b(x)+tau_b(x)); these
# draws permit an apples-to-apples comparison with the Gamma-hurdle benchmark.
hurdle_cate_draws <- pnorm(fit$mu_b + fit$tau_b) - pnorm(fit$mu_b)
hurdle_ate_draws <- rowMeans(hurdle_cate_draws)
hurdle_ate_ci <- quantile(hurdle_ate_draws, c(0.025, 0.975))
saveRDS(list(hurdle_cate_draws = hurdle_cate_draws,
             hurdle_ate_draws = hurdle_ate_draws,
             hurdle_cate_mean = colMeans(hurdle_cate_draws)),
        file.path(outdir, "onc_hurdle_draws.rds"))
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
ate_ci   <- quantile(ate_draws, c(0.025, 0.975))
cat(sprintf("\nATE: %.2f days  95%% CI [%.2f, %.2f]  P(>0)=%.3f\n",
            ate_mean, ate_ci[1], ate_ci[2], mean(ate_draws > 0)))

cate_mean <- colMeans(cate_draws)
res <- data.frame(
  Metric = c("ATE", "Lower 95% CI", "Upper 95% CI", "P(ATE>0)",
             "Mean CATE", "Median CATE"),
  Value  = c(ate_mean, ate_ci[1], ate_ci[2], mean(ate_draws > 0),
             mean(cate_mean), median(cate_mean))
)
write.csv(res, file.path(outdir, "onc_ate_results.csv"), row.names = FALSE)

# ===========================================================================
# Subgroup CATE machinery (identical logic to the MEPS study)
# ===========================================================================
subgroup_cate <- function(idx) {
  sub <- cate_draws[, idx, drop = FALSE]
  draw_means <- rowMeans(sub)
  c(N = ncol(sub),
    CATE = mean(draw_means),
    CI_low  = unname(quantile(draw_means, 0.025)),
    CI_high = unname(quantile(draw_means, 0.975)),
    P_gt_0  = mean(draw_means > 0))
}

contrast_cate <- function(idx_a, idx_b) {
  da <- rowMeans(cate_draws[, idx_a, drop = FALSE])
  db <- rowMeans(cate_draws[, idx_b, drop = FALSE])
  d  <- da - db
  c(Diff = mean(d),
    CI_low  = unname(quantile(d, 0.025)),
    CI_high = unname(quantile(d, 0.975)),
    P_gt_0  = mean(d > 0))
}

# -- interpretable subgroup labels ------------------------------------------
sex_lab    <- factor(df$sex, levels = c("Male", "Female"))
ecog_lab   <- factor(df$b_ecogct, levels = c(0, 1),
                     labels = c("ECOG 0", "ECOG 1"))
diag_lab   <- factor(df$diagtype)
tumcat_lab <- factor(df$tumcat, levels = c("N", "Y"),
                     labels = c("Tumor cat. N", "Tumor cat. Y"))
hpv_lab    <- factor(df$hpv, levels = c("Negative", "Positive", "Unknown"))
dstat_lab  <- factor(df$dstatus,
                     levels = c("newly diagnosed", "recurrent"))
age_lab    <- cut(df$age, breaks = c(-Inf, 49, 59, 69, Inf),
                  labels = c("<50", "50-59", "60-69", "70+"))

group_vars <- list(
  Sex               = sex_lab,
  `ECOG status`     = ecog_lab,
  `Diagnosis site`  = diag_lab,
  `Tumor category`  = tumcat_lab,
  `HPV status`      = hpv_lab,
  `Disease status`  = dstat_lab,
  `Age group`       = age_lab
)

# -- overall + per-subgroup CATE table --------------------------------------
rows <- list()
rows[[1]] <- data.frame(Covariate = "Overall", Level = "All",
                        t(subgroup_cate(rep(TRUE, ncol(cate_draws)))),
                        check.names = FALSE)
for (gname in names(group_vars)) {
  lab <- group_vars[[gname]]
  for (lv in levels(lab)) {
    idx <- which(lab == lv)
    if (length(idx) < 30) next
    rows[[length(rows) + 1]] <- data.frame(
      Covariate = gname, Level = lv,
      t(subgroup_cate(idx)), check.names = FALSE)
  }
}
cate_table <- do.call(rbind, rows)
rownames(cate_table) <- NULL
cate_table[, c("CATE", "CI_low", "CI_high")] <-
  round(cate_table[, c("CATE", "CI_low", "CI_high")], 2)
cate_table$P_gt_0 <- round(cate_table$P_gt_0, 3)
write.csv(cate_table, file.path(outdir, "onc_cate_subgroups.csv"), row.names = FALSE)
cat("\n=== Subgroup CATE estimates ===\n"); print(cate_table)

# -- pairwise heterogeneity contrasts ---------------------------------------
crows <- list()
for (gname in names(group_vars)) {
  lab <- group_vars[[gname]]
  lvs <- levels(lab)[table(lab) >= 30]
  if (length(lvs) < 2) next
  cmb <- combn(lvs, 2)
  for (k in seq_len(ncol(cmb))) {
    a <- cmb[1, k]; b <- cmb[2, k]
    ct <- contrast_cate(which(lab == a), which(lab == b))
    crows[[length(crows) + 1]] <- data.frame(
      Covariate = gname, Contrast = paste(a, "-", b),
      t(ct), check.names = FALSE)
  }
}
contrast_table <- do.call(rbind, crows)
rownames(contrast_table) <- NULL
contrast_table[, c("Diff", "CI_low", "CI_high")] <-
  round(contrast_table[, c("Diff", "CI_low", "CI_high")], 2)
contrast_table$P_gt_0 <- round(contrast_table$P_gt_0, 3)
write.csv(contrast_table, file.path(outdir, "onc_cate_contrasts.csv"), row.names = FALSE)
cat("\n=== Between-group CATE contrasts (heterogeneity) ===\n"); print(contrast_table)

# ===========================================================================
# Figures
# ===========================================================================
XLAB <- "CATE: effect of panitumumab on cumulative grade 3+ AE duration (days)"

# 11 - subgroup CATE forest with 95% credible intervals
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

# 12 - ATE posterior
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

# 13 - unit-level CATE histogram
png(file.path(outdir, "13_ONC_cate_unit_histogram.png"), width = 900, height = 600)
hist(cate_mean, breaks = 50, col = "skyblue", border = "white",
     main = "Distribution of unit-level Conditional Average Treatment Effects",
     xlab = XLAB)
abline(v = ate_mean, col = "firebrick", lwd = 2, lty = 3)
dev.off()

# 14 - subgroup CATE posterior densities (ECOG, tumor category, disease status)
dens_groups <- list(`ECOG status` = ecog_lab,
                    `Tumor category` = tumcat_lab,
                    `Disease status` = dstat_lab)
png(file.path(outdir, "14_ONC_cate_subgroup_densities.png"), width = 1000, height = 850)
par(mfrow = c(2, 2), mar = c(4.5, 4, 3, 1))
for (gname in names(dens_groups)) {
  lab <- dens_groups[[gname]]
  lvs <- levels(lab)[table(lab) >= 30]
  dens <- lapply(lvs, function(lv) density(rowMeans(cate_draws[, which(lab == lv), drop = FALSE])))
  xr <- range(sapply(dens, function(d) range(d$x)))
  yr <- range(sapply(dens, function(d) range(d$y)))
  cols <- c("steelblue", "firebrick", "darkgreen")[seq_along(lvs)]
  plot(NA, xlim = xr, ylim = yr, xlab = "Subgroup CATE (days)", ylab = "Posterior density",
       main = gname)
  for (i in seq_along(dens)) {
    lines(dens[[i]], col = cols[i], lwd = 2)
    polygon(dens[[i]], col = adjustcolor(cols[i], 0.2), border = NA)
  }
  abline(v = 0, lty = 3, col = "grey50")
  legend("topright", legend = lvs, col = cols, lwd = 2, bty = "n")
}
dev.off()

# 15 - unit-level CATE vs age with LOWESS trend
png(file.path(outdir, "15_ONC_cate_vs_age.png"), width = 900, height = 600)
plot(df$age, cate_mean, pch = 19, col = adjustcolor("steelblue", 0.4),
     xlab = "Baseline age (years)", ylab = "Unit posterior-mean CATE (days)",
     main = "Unit-level CATE versus baseline age")
lines(lowess(df$age, cate_mean, f = 0.5), col = "firebrick", lwd = 3)
abline(h = ate_mean, lty = 3, col = "grey40")
dev.off()

# 16 - between-subgroup contrasts forest
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
    " 95% CI [", round(ate_ci[1], 2), ",", round(ate_ci[2], 2), "]\n")
cat("Completed successfully.\n")
