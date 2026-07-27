# Posterior subgroup analysis for the oncology Gamma-hurdle benchmark.
# Run run_gamma_hurdle_oncology.R first.

outdir <- "applied_study_oncology"
draw_file <- file.path(outdir, "gamma_hurdle_draws.rds")
data_file <- file.path(outdir, "zic_bcf_headneck_analysis_data.csv")
if (!file.exists(draw_file)) stop("Missing posterior draws: ", draw_file)
if (!file.exists(data_file)) stop("Missing analysis data: ", data_file)

df <- read.csv(data_file, stringsAsFactors = FALSE)
draws <- readRDS(draw_file)
cate_draws <- draws$cate_draws
hurdle_cate_draws <- draws$hurdle_cate_draws

subgroup_posterior <- function(draws_mat, idx) {
  posterior <- rowMeans(draws_mat[, idx, drop = FALSE])
  c(N = length(idx),
    CATE = mean(posterior),
    CI_low = unname(quantile(posterior, 0.025)),
    CI_high = unname(quantile(posterior, 0.975)),
    P_gt_0 = mean(posterior > 0))
}

contrast_posterior <- function(draws_mat, idx_a, idx_b) {
  diff_draws <- rowMeans(draws_mat[, idx_a, drop = FALSE]) -
    rowMeans(draws_mat[, idx_b, drop = FALSE])
  c(Diff = mean(diff_draws),
    CI_low = unname(quantile(diff_draws, 0.025)),
    CI_high = unname(quantile(diff_draws, 0.975)),
    P_gt_0 = mean(diff_draws > 0))
}

# These labels are intentionally identical to the ZIC-BCF-Smear analysis.
sex_lab <- factor(df$sex, levels = c("Male", "Female"))
ecog_lab <- factor(df$b_ecogct, levels = c(0, 1), labels = c("ECOG 0", "ECOG 1"))
diag_lab <- factor(df$diagtype)
tumcat_lab <- factor(df$tumcat, levels = c("N", "Y"), labels = c("Tumor cat. N", "Tumor cat. Y"))
hpv_lab <- factor(df$hpv, levels = c("Negative", "Positive", "Unknown"))
dstat_lab <- factor(df$dstatus, levels = c("newly diagnosed", "recurrent"))
age_lab <- cut(df$age, breaks = c(-Inf, 49, 59, 69, Inf),
               labels = c("<50", "50-59", "60-69", "70+"))

group_vars <- list(
  Sex = sex_lab,
  `ECOG status` = ecog_lab,
  `Diagnosis site` = diag_lab,
  `Tumor category` = tumcat_lab,
  `HPV status` = hpv_lab,
  `Disease status` = dstat_lab,
  `Age group` = age_lab
)

build_table <- function(draws_mat) {
  all_idx <- seq_len(nrow(df))
  rows <- list(data.frame(Covariate = "Overall", Level = "All",
                          t(subgroup_posterior(draws_mat, all_idx)),
                          check.names = FALSE))
  for (gname in names(group_vars)) {
    lab <- group_vars[[gname]]
    for (level in levels(lab)) {
      idx <- which(lab == level)
      if (length(idx) < 30) next
      rows[[length(rows) + 1]] <- data.frame(
        Covariate = gname, Level = level,
        t(subgroup_posterior(draws_mat, idx)), check.names = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[, c("CATE", "CI_low", "CI_high")] <- round(out[, c("CATE", "CI_low", "CI_high")], 4)
  out$P_gt_0 <- round(out$P_gt_0, 3)
  out
}

build_contrasts <- function(draws_mat) {
  rows <- list()
  for (gname in names(group_vars)) {
    lab <- group_vars[[gname]]
    eligible <- levels(lab)[table(lab) >= 30]
    if (length(eligible) < 2) next
    pairs <- combn(eligible, 2)
    for (k in seq_len(ncol(pairs))) {
      a <- pairs[1, k]
      b <- pairs[2, k]
      contrast <- contrast_posterior(draws_mat, which(lab == a), which(lab == b))
      rows[[length(rows) + 1]] <- data.frame(
        Covariate = gname, Contrast = paste(a, "-", b),
        t(contrast), check.names = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[, c("Diff", "CI_low", "CI_high")] <- round(out[, c("Diff", "CI_low", "CI_high")], 4)
  out$P_gt_0 <- round(out$P_gt_0, 3)
  out
}

cate_table <- build_table(cate_draws)
hurdle_table <- build_table(hurdle_cate_draws)
contrast_table <- build_contrasts(cate_draws)
hurdle_contrast_table <- build_contrasts(hurdle_cate_draws)

write.csv(cate_table, file.path(outdir, "gamma_hurdle_cate_subgroups.csv"), row.names = FALSE)
write.csv(hurdle_table, file.path(outdir, "gamma_hurdle_hurdle_subgroups.csv"), row.names = FALSE)
write.csv(contrast_table, file.path(outdir, "gamma_hurdle_cate_contrasts.csv"), row.names = FALSE)
write.csv(hurdle_contrast_table, file.path(outdir, "gamma_hurdle_hurdle_contrasts.csv"), row.names = FALSE)

plot_df <- cate_table[cate_table$Covariate != "Overall", ]
plot_df$row_lab <- paste0(plot_df$Covariate, ": ", plot_df$Level)
plot_df <- plot_df[nrow(plot_df):1, ]
overall <- cate_table$CATE[cate_table$Covariate == "Overall"]

png(file.path(outdir, "17_ONC_gamma_hurdle_cate_subgroups.png"), width = 1000, height = 1100)
par(mar = c(5, 16, 4, 2))
yy <- seq_len(nrow(plot_df))
plot(plot_df$CATE, yy, xlim = range(c(plot_df$CI_low, plot_df$CI_high, 0)),
     yaxt = "n", pch = 19, col = "steelblue", cex = 1.3,
     xlab = "CATE: effect of panitumumab on cumulative grade 3+ AE duration (days)",
     ylab = "", main = "Gamma Hurdle: subgroup CATE estimates with 95% CrI")
axis(2, at = yy, labels = plot_df$row_lab, las = 1, cex.axis = 0.9)
segments(plot_df$CI_low, yy, plot_df$CI_high, yy, col = "steelblue", lwd = 2)
abline(v = 0, lty = 2, col = "grey50")
abline(v = overall, lty = 3, col = "firebrick", lwd = 2)
legend("bottomright", legend = c("Subgroup CATE (95% CrI)", "No effect",
                                sprintf("Overall ATE = %.1f days", overall)),
       col = c("steelblue", "grey50", "firebrick"), pch = c(19, NA, NA),
       lty = c(NA, 2, 3), lwd = c(2, 1, 2), bty = "n")
dev.off()

cat("Completed Gamma-hurdle posterior subgroup analysis.\n")
