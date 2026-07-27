# Create forest plots supporting the ZIC-BCF distinct-calendar-day CATE result.
# This script uses saved posterior summaries, so it does not refit either model.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

script_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))])
if (length(script_file) != 1L) stop("Run this script with Rscript.")

data_dir <- dirname(normalizePath(script_file))
out_dir <- file.path(data_dir, "analysis")

subgroups <- read_csv(file.path(out_dir, "nonoverlap_subgroup_model_results.csv"), show_col_types = FALSE) %>%
  filter(
    Model == "ZIC-BCF-Smear",
    Estimand == "Total non-overlapping severe-AE days",
    Covariate != "Overall"
  ) %>%
  mutate(label = paste(Covariate, Level, sep = ": "))

contrasts <- read_csv(file.path(out_dir, "nonoverlap_subgroup_contrasts.csv"), show_col_types = FALSE) %>%
  filter(
    Model == "ZIC-BCF-Smear",
    Estimand == "Total non-overlapping severe-AE days"
  ) %>%
  mutate(
    label = paste(Covariate, Contrast, sep = ": "),
    credible = CI_low > 0 | CI_high < 0
  )

if (nrow(subgroups) == 0L || nrow(contrasts) == 0L) {
  stop("Missing ZIC-BCF subgroup results. Run run_nonoverlap_calendar_day_analysis.R first.")
}

# Overall ATE reference line, matching the summed-duration subgroup figures.
overall_ate <- read_csv(file.path(out_dir, "nonoverlap_overall_model_results.csv"), show_col_types = FALSE) %>%
  filter(Model == "ZIC-BCF-Smear", Estimand == "Total non-overlapping severe-AE days") %>%
  pull(Estimate)

if (length(overall_ate) != 1L) stop("Could not read a unique ZIC-BCF-Smear overall ATE.")

# Figure 1: treatment effects within baseline subgroups.  These intervals
# answer whether the treatment effect is positive in each subgroup, whereas
# the contrast plot below answers whether effects differ between subgroups.
subgroups <- subgroups[nrow(subgroups):1L, ]
png(file.path(out_dir, "nonoverlap_zic_cate_subgroups.png"), width = 1150, height = 1150)
par(mar = c(5, 20, 4, 2))
positions <- seq_len(nrow(subgroups))
plot(
  subgroups$Estimate, positions,
  xlim = range(c(subgroups$CI_low, subgroups$CI_high, 0)),
  ylim = c(0.5, nrow(subgroups) + 0.5), yaxt = "n", pch = 19,
  col = "steelblue", cex = 1.1,
  xlab = "CATE: effect of panitumumab on distinct severe-AE days",
  ylab = "", main = "ZIC-BCF subgroup CATEs with 95% credible intervals"
)
segments(subgroups$CI_low, positions, subgroups$CI_high, positions,
         col = "steelblue", lwd = 2)
axis(2, at = positions, labels = subgroups$label, las = 1, cex.axis = 0.85)
abline(v = 0, lty = 2, col = "grey40")
abline(v = overall_ate, lty = 3, col = "darkred")
legend("bottomright",
       legend = c("Subgroup CATE (95% CrI)", "No effect",
                  sprintf("Overall ATE = %.1f days", overall_ate)),
       col = c("steelblue", "grey40", "darkred"), pch = c(19, NA, NA),
       lty = c(NA, 2, 3), lwd = c(2, 1, 1), bty = "n")
dev.off()

# Figure 2: direct evidence for heterogeneity.  Every interval crosses zero,
# so no prespecified contrast is credible on the 95% CrI criterion.
contrasts <- contrasts[nrow(contrasts):1L, ]
png(file.path(out_dir, "nonoverlap_zic_cate_contrasts.png"), width = 1200, height = 1100)
par(mar = c(5, 24, 4, 2))
positions <- seq_len(nrow(contrasts))
point_colours <- ifelse(contrasts$credible, "firebrick", "grey45")
plot(
  contrasts$Estimate, positions,
  xlim = range(c(contrasts$CI_low, contrasts$CI_high, 0)),
  ylim = c(0.5, nrow(contrasts) + 0.5), yaxt = "n", pch = 19,
  col = point_colours, cex = 1.1,
  xlab = "Between-subgroup difference in CATE (days)", ylab = "",
  main = "ZIC-BCF subgroup CATE contrasts with 95% credible intervals"
)
segments(contrasts$CI_low, positions, contrasts$CI_high, positions,
         col = point_colours, lwd = 2)
axis(2, at = positions, labels = contrasts$label, las = 1, cex.axis = 0.78)
abline(v = 0, lty = 2, col = "grey40")
legend("bottomright",
       legend = c("Credible contrast (CrI excludes 0)", "Contrast not credible"),
       col = c("firebrick", "grey45"), pch = 19, bty = "n")
dev.off()

cat("Created ZIC-BCF subgroup and contrast forest plots.\n")
cat("Credible total-outcome contrasts:", sum(contrasts$credible), "of", nrow(contrasts), "\n")
