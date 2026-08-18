# Subgroup analysis for the Gamma-hurdle benchmark on the MEPS dental sample.
#
# Run run_gamma_hurdle.R first; this script consumes the posterior summaries it
# saves and performs no fitting. Subgroup effects, contrasts and the two
# estimands (analytic-sample average and survey-weighted population target) are
# constructed exactly as in run_cate_subgroups.R, so the two models are always
# compared on identical quantities.

source("applied_study/meps_common.R")

OUTDIR <- "applied_study"
draw_file <- file.path(OUTDIR, "gamma_hurdle_draws.rds")
if (!file.exists(draw_file)) stop("Run run_gamma_hurdle.R first.")

posterior <- readRDS(draw_file)
df <- meps_load(file.path(OUTDIR, "h251.csv"))
labels <- meps_subgroup_labels(df)
weights <- df$PERWT23F

domain_names <- posterior$domain_names
covariate <- vapply(posterior$domains, function(d) d$covariate, character(1))
level <- vapply(posterior$domains, function(d) d$level, character(1))

build_subgroup_table <- function(draws, design_variance, scale_label) {
  do.call(rbind, lapply(seq_along(domain_names), function(k) {
    unweighted <- meps_posterior_summary(draws$unweighted[, k])
    weighted <- meps_posterior_summary(draws$weighted[, k], design_variance[k])
    data.frame(
      Scale = scale_label, Covariate = covariate[k], Level = level[k],
      N = posterior$domain_n[k],
      Sample_estimate = unweighted[["Estimate"]],
      Sample_CI_low = unweighted[["CI_low"]],
      Sample_CI_high = unweighted[["CI_high"]],
      Population_estimate = weighted[["Estimate"]],
      Population_CI_low_design = weighted[["CI_low_design"]],
      Population_CI_high_design = weighted[["CI_high_design"]],
      check.names = FALSE
    )
  }))
}

build_contrast_table <- function(draws, unit_mean, scale_label) {
  rows <- list()
  for (group_name in setdiff(unique(covariate), "Overall")) {
    members <- which(covariate == group_name)
    if (length(members) < 2L) next
    comparisons <- combn(members, 2L)
    for (column in seq_len(ncol(comparisons))) {
      a <- comparisons[1L, column]
      b <- comparisons[2L, column]
      label <- labels[[group_name]]
      design_variance <- meps_design_variance_contrast(
        unit_mean, weights, df$VARSTR, df$VARPSU,
        !is.na(label) & label == level[a], !is.na(label) & label == level[b])
      sample_summary <- meps_posterior_summary(draws$unweighted[, a] - draws$unweighted[, b])
      population_summary <- meps_posterior_summary(
        draws$weighted[, a] - draws$weighted[, b], design_variance)
      rows[[length(rows) + 1L]] <- data.frame(
        Scale = scale_label, Covariate = group_name,
        Contrast = paste(level[a], "-", level[b]),
        Sample_estimate = sample_summary[["Estimate"]],
        Sample_CI_low = sample_summary[["CI_low"]],
        Sample_CI_high = sample_summary[["CI_high"]],
        Population_estimate = population_summary[["Estimate"]],
        Population_CI_low_design = population_summary[["CI_low_design"]],
        Population_CI_high_design = population_summary[["CI_high_design"]],
        check.names = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

cate_table <- build_subgroup_table(posterior$cate, posterior$cate_design_variance,
                                   "Dollar-scale CATE")
hurdle_table <- build_subgroup_table(posterior$hurdle, posterior$hurdle_design_variance,
                                     "Participation-margin CATE")
cate_contrasts <- build_contrast_table(posterior$cate, posterior$cate_unit_mean,
                                       "Dollar-scale CATE")
hurdle_contrasts <- build_contrast_table(posterior$hurdle, posterior$hurdle_unit_mean,
                                         "Participation-margin CATE")

write.csv(cate_table, file.path(OUTDIR, "gamma_hurdle_cate_subgroups.csv"), row.names = FALSE)
write.csv(hurdle_table, file.path(OUTDIR, "gamma_hurdle_hurdle_subgroups.csv"), row.names = FALSE)
write.csv(cate_contrasts, file.path(OUTDIR, "gamma_hurdle_cate_contrasts.csv"), row.names = FALSE)
write.csv(hurdle_contrasts, file.path(OUTDIR, "gamma_hurdle_hurdle_contrasts.csv"), row.names = FALSE)

plot_df <- cate_table[cate_table$Covariate != "Overall", ]
plot_df$row_lab <- paste0(plot_df$Covariate, ": ", plot_df$Level)
plot_df <- plot_df[nrow(plot_df):1, ]
overall <- cate_table$Sample_estimate[cate_table$Covariate == "Overall"]

png(file.path(OUTDIR, "gamma_hurdle_cate_subgroups.png"), width = 1000, height = 1100)
par(mar = c(5, 20, 4, 2))
yy <- seq_len(nrow(plot_df))
plot(plot_df$Sample_estimate, yy,
     xlim = range(c(plot_df$Sample_CI_low, plot_df$Sample_CI_high, 0)),
     yaxt = "n", pch = 19, col = "steelblue", cex = 1.3,
     xlab = "CATE: effect of dental insurance on annual dental expenditure ($)",
     ylab = "", main = "Gamma Hurdle: subgroup CATEs with 95% credible intervals")
axis(2, at = yy, labels = plot_df$row_lab, las = 1, cex.axis = 0.85)
segments(plot_df$Sample_CI_low, yy, plot_df$Sample_CI_high, yy, col = "steelblue", lwd = 2)
abline(v = 0, lty = 2, col = "grey50")
abline(v = overall, lty = 3, col = "firebrick", lwd = 2)
legend("bottomright",
       legend = c("Subgroup CATE (95% CrI)", "No effect",
                  sprintf("Overall ATE = $%.0f", overall)),
       col = c("steelblue", "grey50", "firebrick"), pch = c(19, NA, NA),
       lty = c(NA, 2, 3), lwd = c(2, 1, 2), bty = "n")
dev.off()

cat("=== Gamma hurdle: dollar-scale subgroup CATEs ===\n")
print(cate_table[, c("Covariate", "Level", "N", "Sample_estimate", "Sample_CI_low",
                     "Sample_CI_high", "Population_estimate")],
      row.names = FALSE, digits = 4)
cat("\nCompleted Gamma-hurdle posterior subgroup analysis.\n")
