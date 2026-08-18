# ==========================================================================
# Subgroup (CATE) analysis for the MEPS dental-insurance study.
#
# Run run_study.R first; this script consumes the posterior summaries it saves
# and performs no fitting.
#
# A subgroup CATE is the average of unit-level effects over the units in the
# subgroup, taken WITHIN each posterior draw, which yields a posterior
# distribution for the subgroup effect and therefore a point estimate and a 95%
# credible interval. Between-group contrasts are formed the same way, so a
# contrast inherits the posterior covariance of its two subgroups rather than
# being read off two unrelated marginal intervals.
#
# TWO ESTIMANDS ARE REPORTED FOR EVERY QUANTITY.
#
#   Analytic-sample average. The unweighted average over the 18,763 records.
#   This is a sample quantity over an oversampled design, not a United States
#   population quantity.
#
#   Survey-weighted population target. The average weighted by PERWT23F, with a
#   design-based variance component from the stratification (VARSTR) and
#   clustering (VARPSU) added to the posterior variance. Only this column
#   supports a national reading.
#
# The two are reported side by side because they answer different questions and
# because the difference between them is itself informative about how far the
# oversampled subgroups drive the sample-average result.
# ==========================================================================

source("applied_study/meps_common.R")

OUTDIR <- "applied_study"
summary_file <- file.path(OUTDIR, "meps_zicbcf_posterior_summaries.rds")
if (!file.exists(summary_file)) stop("Run run_study.R first.")

posterior <- readRDS(summary_file)
df <- meps_load(file.path(OUTDIR, "h251.csv"))
labels <- meps_subgroup_labels(df)
weights <- df$PERWT23F

domain_names <- posterior$domain_names
covariate <- vapply(posterior$domains, function(d) d$covariate, character(1))
level <- vapply(posterior$domains, function(d) d$level, character(1))

cat(sprintf("Posterior: %d chains, %d retained draws each, %d pooled draws.\n",
            posterior$n_chains, posterior$n_sim, posterior$total_draws))

# ---------------------------------------------------------------------------
# Subgroup tables
# ---------------------------------------------------------------------------
build_subgroup_table <- function(draws, design_variance, scale_label) {
  rows <- lapply(seq_along(domain_names), function(k) {
    unweighted <- meps_posterior_summary(draws$unweighted[, k])
    weighted <- meps_posterior_summary(draws$weighted[, k], design_variance[k])
    data.frame(
      Scale = scale_label,
      Covariate = covariate[k],
      Level = level[k],
      N = posterior$domain_n[k],
      Weighted_N = posterior$domain_weighted_n[k],
      Sample_estimate = unweighted[["Estimate"]],
      Sample_CI_low = unweighted[["CI_low"]],
      Sample_CI_high = unweighted[["CI_high"]],
      Population_estimate = weighted[["Estimate"]],
      Population_CI_low = weighted[["CI_low"]],
      Population_CI_high = weighted[["CI_high"]],
      Population_design_SE = weighted[["Design_SE"]],
      Population_total_SE = weighted[["Total_SE"]],
      Population_CI_low_design = weighted[["CI_low_design"]],
      Population_CI_high_design = weighted[["CI_high_design"]],
      P_gt_0 = weighted[["P_gt_0"]],
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

cate_table <- build_subgroup_table(posterior$cate, posterior$cate_design_variance,
                                   "Dollar-scale CATE")
hurdle_table <- build_subgroup_table(posterior$hurdle, posterior$hurdle_design_variance,
                                     "Participation-margin CATE")

write.csv(cate_table, file.path(OUTDIR, "cate_subgroups.csv"), row.names = FALSE)
write.csv(hurdle_table, file.path(OUTDIR, "zicbcf_hurdle_subgroups.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Pairwise contrasts within each covariate
# ---------------------------------------------------------------------------
build_contrast_table <- function(draws, unit_mean, scale_label) {
  rows <- list()
  for (group_name in setdiff(unique(covariate), "Overall")) {
    members <- which(covariate == group_name)
    if (length(members) < 2L) next
    comparisons <- combn(members, 2L)
    for (column in seq_len(ncol(comparisons))) {
      a <- comparisons[1L, column]
      b <- comparisons[2L, column]
      unweighted <- draws$unweighted[, a] - draws$unweighted[, b]
      weighted <- draws$weighted[, a] - draws$weighted[, b]

      label <- labels[[group_name]]
      domain_a <- !is.na(label) & label == level[a]
      domain_b <- !is.na(label) & label == level[b]
      design_variance <- meps_design_variance_contrast(
        unit_mean, weights, df$VARSTR, df$VARPSU, domain_a, domain_b)

      sample_summary <- meps_posterior_summary(unweighted)
      population_summary <- meps_posterior_summary(weighted, design_variance)
      rows[[length(rows) + 1L]] <- data.frame(
        Scale = scale_label,
        Covariate = group_name,
        Contrast = paste(level[a], "-", level[b]),
        Sample_estimate = sample_summary[["Estimate"]],
        Sample_CI_low = sample_summary[["CI_low"]],
        Sample_CI_high = sample_summary[["CI_high"]],
        Population_estimate = population_summary[["Estimate"]],
        Population_CI_low = population_summary[["CI_low"]],
        Population_CI_high = population_summary[["CI_high"]],
        Population_design_SE = population_summary[["Design_SE"]],
        Population_CI_low_design = population_summary[["CI_low_design"]],
        Population_CI_high_design = population_summary[["CI_high_design"]],
        P_gt_0 = population_summary[["P_gt_0"]],
        check.names = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

cate_contrasts <- build_contrast_table(posterior$cate, posterior$cate_unit_mean,
                                       "Dollar-scale CATE")
hurdle_contrasts <- build_contrast_table(posterior$hurdle, posterior$hurdle_unit_mean,
                                         "Participation-margin CATE")

write.csv(cate_contrasts, file.path(OUTDIR, "cate_contrasts.csv"), row.names = FALSE)
write.csv(hurdle_contrasts, file.path(OUTDIR, "zicbcf_hurdle_contrasts.csv"),
          row.names = FALSE)

credible <- function(table, low, high) {
  table[table[[low]] > 0 | table[[high]] < 0, ]
}
cat("\n=== Dollar-scale contrasts credible on the analytic-sample average ===\n")
print(credible(cate_contrasts, "Sample_CI_low", "Sample_CI_high")[
  , c("Covariate", "Contrast", "Sample_estimate", "Sample_CI_low", "Sample_CI_high")],
  row.names = FALSE, digits = 4)
cat("\n=== Dollar-scale contrasts credible on the design-aware population target ===\n")
print(credible(cate_contrasts, "Population_CI_low_design", "Population_CI_high_design")[
  , c("Covariate", "Contrast", "Population_estimate", "Population_CI_low_design",
      "Population_CI_high_design")], row.names = FALSE, digits = 4)

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
XLAB <- "CATE: effect of dental insurance on annual dental expenditure ($)"
overall <- which(domain_names == "Overall | All")

plot_forest <- function(table, file, estimate, low, high, xlab, title, overall_value) {
  plot_df <- table[table$Covariate != "Overall", ]
  plot_df$row_lab <- paste0(plot_df$Covariate, ": ", plot_df$Level)
  plot_df <- plot_df[nrow(plot_df):1, ]
  png(file, width = 1000, height = 1100)
  par(mar = c(5, 20, 4, 2))
  yy <- seq_len(nrow(plot_df))
  plot(plot_df[[estimate]], yy,
       xlim = range(c(plot_df[[low]], plot_df[[high]], 0)),
       yaxt = "n", pch = 19, col = "steelblue", cex = 1.3,
       xlab = xlab, ylab = "", main = title)
  axis(2, at = yy, labels = plot_df$row_lab, las = 1, cex.axis = 0.85)
  segments(plot_df[[low]], yy, plot_df[[high]], yy, col = "steelblue", lwd = 2)
  abline(v = 0, lty = 2, col = "grey50")
  abline(v = overall_value, lty = 3, col = "firebrick", lwd = 2)
  legend("bottomright",
         legend = c("Subgroup CATE (95% interval)", "No effect",
                    sprintf("Overall = %.0f", overall_value)),
         col = c("steelblue", "grey50", "firebrick"), pch = c(19, NA, NA),
         lty = c(NA, 2, 3), lwd = c(2, 1, 2), bty = "n")
  dev.off()
}

plot_forest(cate_table, file.path(OUTDIR, "11_MEPS_cate_subgroups.png"),
            "Sample_estimate", "Sample_CI_low", "Sample_CI_high", XLAB,
            "Subgroup CATEs, analytic-sample average, with 95% credible intervals",
            cate_table$Sample_estimate[overall])
plot_forest(cate_table, file.path(OUTDIR, "11b_MEPS_cate_subgroups_weighted.png"),
            "Population_estimate", "Population_CI_low_design", "Population_CI_high_design",
            XLAB,
            "Subgroup CATEs, survey-weighted population target, design-aware intervals",
            cate_table$Population_estimate[overall])

ate_draws <- posterior$cate$unweighted[, overall]
weighted_ate_draws <- posterior$cate$weighted[, overall]
png(file.path(OUTDIR, "12_MEPS_ate_posterior.png"), width = 900, height = 600)
hist(weighted_ate_draws, breaks = 40, col = "skyblue", border = "white",
     main = "Posterior of the survey-weighted average treatment effect",
     xlab = XLAB)
abline(v = mean(weighted_ate_draws), col = "firebrick", lwd = 2)
abline(v = c(cate_table$Population_CI_low_design[overall],
             cate_table$Population_CI_high_design[overall]),
       col = "firebrick", lwd = 2, lty = 2)
legend("topright",
       legend = c(sprintf("Posterior mean = $%.0f", mean(weighted_ate_draws)),
                  sprintf("Design-aware 95%% interval [$%.0f, $%.0f]",
                          cate_table$Population_CI_low_design[overall],
                          cate_table$Population_CI_high_design[overall])),
       col = "firebrick", lwd = 2, lty = c(1, 2), bty = "n")
dev.off()

unit_level <- readRDS(file.path(OUTDIR, "meps_unit_level_posterior_means.rds"))
png(file.path(OUTDIR, "13_MEPS_cate_unit_histogram.png"), width = 900, height = 600)
hist(unit_level$cate_posterior_mean, breaks = 50, col = "skyblue", border = "white",
     main = "Distribution of unit-level Conditional Average Treatment Effects",
     xlab = XLAB)
abline(v = mean(ate_draws), col = "firebrick", lwd = 2, lty = 3)
dev.off()

density_groups <- c("Race and ethnicity", "Age group")
png(file.path(OUTDIR, "14_MEPS_cate_subgroup_densities.png"), width = 1000, height = 500)
par(mfrow = c(1, 2), mar = c(4.5, 4, 3, 1))
for (group_name in density_groups) {
  members <- which(covariate == group_name)
  densities <- lapply(members, function(k) density(posterior$cate$weighted[, k]))
  xr <- range(sapply(densities, function(d) range(d$x)))
  yr <- range(sapply(densities, function(d) range(d$y)))
  cols <- hcl.colors(length(members), "Dark 3")
  plot(NA, xlim = xr, ylim = yr, xlab = "Subgroup CATE ($)",
       ylab = "Posterior density", main = group_name)
  for (i in seq_along(densities)) {
    lines(densities[[i]], col = cols[i], lwd = 2)
    polygon(densities[[i]], col = adjustcolor(cols[i], 0.2), border = NA)
  }
  legend("topright", legend = level[members], col = cols, lwd = 2, bty = "n",
         cex = 0.8)
}
dev.off()

png(file.path(OUTDIR, "15_MEPS_cate_vs_covariates.png"), width = 1000, height = 500)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
plot(unit_level$AGE23X, unit_level$cate_posterior_mean, pch = 19,
     col = adjustcolor("steelblue", 0.2), xlab = "Age (years)",
     ylab = "Unit posterior-mean CATE ($)", main = "CATE versus age")
lines(lowess(unit_level$AGE23X, unit_level$cate_posterior_mean, f = 0.4),
      col = "firebrick", lwd = 3)
plot(unit_level$FAMINC23, unit_level$cate_posterior_mean, pch = 19,
     col = adjustcolor("steelblue", 0.2), xlab = "Family income ($)",
     ylab = "Unit posterior-mean CATE ($)", main = "CATE versus family income")
lines(lowess(unit_level$FAMINC23, unit_level$cate_posterior_mean, f = 0.4),
      col = "firebrick", lwd = 3)
dev.off()

contrast_plot <- cate_contrasts
contrast_plot$row_lab <- paste0(contrast_plot$Covariate, ": ", contrast_plot$Contrast)
contrast_plot <- contrast_plot[nrow(contrast_plot):1, ]
is_credible <- contrast_plot$Population_CI_low_design > 0 |
  contrast_plot$Population_CI_high_design < 0
png(file.path(OUTDIR, "16_MEPS_cate_contrasts_forest.png"), width = 1000, height = 1300)
par(mar = c(5, 22, 4, 2))
yy <- seq_len(nrow(contrast_plot))
cols <- ifelse(is_credible, "firebrick", "grey50")
plot(contrast_plot$Population_estimate, yy,
     xlim = range(c(contrast_plot$Population_CI_low_design,
                    contrast_plot$Population_CI_high_design, 0)),
     yaxt = "n", pch = 19, col = cols, cex = 1.2,
     xlab = "Between-subgroup difference in CATE ($), survey-weighted",
     ylab = "", main = "Pairwise subgroup contrasts with design-aware 95% intervals")
axis(2, at = yy, labels = contrast_plot$row_lab, las = 1, cex.axis = 0.7)
segments(contrast_plot$Population_CI_low_design, yy,
         contrast_plot$Population_CI_high_design, yy, col = cols, lwd = 2)
abline(v = 0, lty = 2, col = "grey40")
legend("bottomright", legend = c("Credible (interval excludes 0)", "Not credible"),
       col = c("firebrick", "grey50"), pch = 19, bty = "n")
dev.off()

cat("\nCompleted successfully.\n")
