# ==========================================================================
# Regenerate the oncology subgroup tables, contrast tables and subgroup
# figures from the saved posterior draws, after dropping the composite
# prior-treatment indicator PRHNTRTC ("Prior SCCHN treatment") from the
# prespecified subgroup list.
#
# Review point: PRRADIO and PRHNTRTC differ for four patients, so every
# subgroup grid previously carried a near-duplicate contrast per endpoint per
# model in an analysis that is already unadjusted for multiplicity. The
# subgroup list is now defined once in oncology_common.R without the
# composite; this script rebuilds all downstream artifacts from the committed
# posterior draws so the numbers in the tables are unchanged except for the
# removal of the duplicated strata. No refitting takes place.
#
# Regenerated:
#   onc_cate_subgroups.csv, onc_cate_contrasts.csv            (summed, ZIC)
#   11_ONC_cate_subgroups.png, 14_ONC_cate_subgroup_densities.png
#   gamma_hurdle_{cate,hurdle}_{subgroups,contrasts}.csv      (summed, Gamma)
#   17_ONC_gamma_hurdle_cate_subgroups.png
#   New_Data/analysis/nonoverlap_subgroup_model_results.csv,
#   New_Data/analysis/nonoverlap_subgroup_contrasts.csv       (distinct days)
# ==========================================================================

outdir <- "applied_study_oncology"
source(file.path(outdir, "oncology_common.R"))

data_file <- file.path(outdir, "zic_bcf_headneck_analysis_data.csv")
if (!file.exists(data_file)) stop("Run prepare_analysis_data.R first.")
df <- read.csv(data_file, stringsAsFactors = FALSE)
stopifnot(nrow(df) == 520L)

posterior_summary <- function(draws) {
  ci <- stats::quantile(draws, c(0.025, 0.975))
  c(Estimate = mean(draws), CI_low = unname(ci[1L]), CI_high = unname(ci[2L]),
    P_gt_0 = mean(draws > 0))
}

build_table <- function(draws_mat, model, estimand) {
  rows <- list(data.frame(Model = model, Estimand = estimand,
                          Covariate = "Overall", Level = "All",
                          N = ncol(draws_mat),
                          t(posterior_summary(rowMeans(draws_mat))),
                          check.names = FALSE))
  labels <- onc_subgroup_labels(df)
  for (group_name in names(labels)) {
    lab <- labels[[group_name]]
    for (level in levels(lab)) {
      idx <- which(lab == level)
      if (length(idx) < ONC_MIN_SUBGROUP_N) next
      rows[[length(rows) + 1L]] <- data.frame(
        Model = model, Estimand = estimand, Covariate = group_name,
        Level = level, N = length(idx),
        t(posterior_summary(rowMeans(draws_mat[, idx, drop = FALSE]))),
        check.names = FALSE)
    }
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

build_contrasts <- function(draws_mat, model, estimand) {
  labels <- onc_subgroup_labels(df)
  rows <- list()
  for (group_name in names(labels)) {
    lab <- labels[[group_name]]
    eligible <- levels(lab)[table(lab) >= ONC_MIN_SUBGROUP_N]
    if (length(eligible) < 2L) next
    pairs <- combn(eligible, 2L)
    for (k in seq_len(ncol(pairs))) {
      a <- pairs[1, k]; b <- pairs[2, k]
      diff_draws <- rowMeans(draws_mat[, which(lab == a), drop = FALSE]) -
        rowMeans(draws_mat[, which(lab == b), drop = FALSE])
      rows[[length(rows) + 1L]] <- data.frame(
        Model = model, Estimand = estimand, Covariate = group_name,
        Contrast = paste(a, "-", b), t(posterior_summary(diff_draws)),
        check.names = FALSE)
    }
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

cate_forest_plot <- function(table, file, title) {
  names(table)[names(table) == "CATE"] <- "Estimate"
  plot_df <- table[table$Covariate != "Overall", ]
  plot_df$row_lab <- paste0(plot_df$Covariate, ": ", plot_df$Level)
  plot_df <- plot_df[nrow(plot_df):1, ]
  overall <- table$Estimate[table$Covariate == "Overall"]
  grDevices::png(file.path(outdir, file), width = 1000, height = 1100)
  par(mar = c(5, 16, 4, 2))
  yy <- seq_len(nrow(plot_df))
  plot(plot_df$Estimate, yy, xlim = range(c(plot_df$CI_low, plot_df$CI_high, 0)),
       yaxt = "n", pch = 19, col = "steelblue", cex = 1.3, ylab = "",
       xlab = "CATE: effect of panitumumab on cumulative grade 3+ AE duration (days)",
       main = title)
  graphics::axis(2, at = yy, labels = plot_df$row_lab, las = 1, cex.axis = 0.9)
  graphics::segments(plot_df$CI_low, yy, plot_df$CI_high, yy,
                     col = "steelblue", lwd = 2)
  graphics::abline(v = 0, lty = 2, col = "grey50")
  graphics::abline(v = overall, lty = 3, col = "firebrick", lwd = 2)
  graphics::legend("bottomright",
                   legend = c("Subgroup CATE (95% CrI)", "No effect",
                              sprintf("Overall ATE = %.1f days", overall)),
                   col = c("steelblue", "grey50", "firebrick"),
                   pch = c(19, NA, NA), lty = c(NA, 2, 3),
                   lwd = c(2, 1, 2), bty = "n")
  grDevices::dev.off()
}

# ---- 1. Summed duration, ZIC-BCF-Smear ------------------------------------
cat("ZIC-BCF-Smear, summed duration...\n")
cate_draws <- readRDS(file.path(outdir, "onc_cate_draws.rds"))
stopifnot(ncol(cate_draws) == nrow(df))
zic_table <- build_table(cate_draws, "ZIC-BCF-Smear",
                         "Summed severe-AE-record duration (event-days)")
names(zic_table)[names(zic_table) == "Estimate"] <- "CATE"
zic_table[, c("CATE", "CI_low", "CI_high")] <-
  round(zic_table[, c("CATE", "CI_low", "CI_high")], 2)
zic_table$P_gt_0 <- round(zic_table$P_gt_0, 3)
write.csv(zic_table, file.path(outdir, "onc_cate_subgroups.csv"), row.names = FALSE)
zc <- build_contrasts(cate_draws, "ZIC-BCF-Smear",
                      "Summed severe-AE-record duration (event-days)")
names(zc)[names(zc) == "Estimate"] <- "Diff"
zc[, c("Diff", "CI_low", "CI_high")] <-
  round(zc[, c("Diff", "CI_low", "CI_high")], 2)
zc$P_gt_0 <- round(zc$P_gt_0, 3)
write.csv(zc, file.path(outdir, "onc_cate_contrasts.csv"), row.names = FALSE)
cate_forest_plot(zic_table, "11_ONC_cate_subgroups.png",
                 "Subgroup CATE estimates with 95% credible intervals")

labels <- onc_subgroup_labels(df)
dens_groups <- labels[c("ECOG status", "Disease status")]
png(file.path(outdir, "14_ONC_cate_subgroup_densities.png"), width = 1000, height = 500)
par(mfrow = c(1, 2), mar = c(4.5, 4, 3, 1))
for (gname in names(dens_groups)) {
  lab <- dens_groups[[gname]]
  lvs <- levels(lab)[table(lab) >= ONC_MIN_SUBGROUP_N]
  dens <- lapply(lvs, function(lv)
    density(rowMeans(cate_draws[, which(lab == lv), drop = FALSE])))
  xr <- range(sapply(dens, function(d) range(d$x)))
  plot(NA, xlim = xr, ylim = range(sapply(dens, function(d) max(d$y))),
       main = gname, xlab = "CATE (days)")
  for (i in seq_along(lvs)) lines(dens[[i]], col = i + 1, lwd = 2)
  legend("topright", legend = lvs, col = seq_along(lvs) + 1, lwd = 2, bty = "n")
}
dev.off()

# ---- 2. Summed duration, Gamma-hurdle benchmark -----------------------------
cat("Gamma hurdle, summed duration...\n")
gamma <- readRDS(file.path(outdir, "gamma_hurdle_draws.rds"))
stopifnot(ncol(gamma$cate_draws) == nrow(df))
g_round4 <- function(t) {
  t[, c("Estimate", "CI_low", "CI_high")] <-
    round(t[, c("Estimate", "CI_low", "CI_high")], 4)
  t$P_gt_0 <- round(t$P_gt_0, 3); t
}
gt <- build_table(gamma$cate_draws, "", "")
gt <- g_round4(gt); names(gt)[names(gt) == "Estimate"] <- "CATE"
write.csv(gt[, c("Covariate","Level","N","CATE","CI_low","CI_high","P_gt_0")],
          file.path(outdir, "gamma_hurdle_cate_subgroups.csv"), row.names = FALSE)
ht <- build_table(gamma$hurdle_cate_draws, "", "")
ht <- g_round4(ht); names(ht)[names(ht) == "Estimate"] <- "CATE"
write.csv(ht[, c("Covariate","Level","N","CATE","CI_low","CI_high","P_gt_0")],
          file.path(outdir, "gamma_hurdle_hurdle_subgroups.csv"), row.names = FALSE)
gc_t <- build_contrasts(gamma$cate_draws, "", "")
gc_t <- g_round4(gc_t); names(gc_t)[names(gc_t) == "Estimate"] <- "Diff"
write.csv(gc_t[, c("Covariate","Contrast","Diff","CI_low","CI_high","P_gt_0")],
          file.path(outdir, "gamma_hurdle_cate_contrasts.csv"), row.names = FALSE)
gh_t <- build_contrasts(gamma$hurdle_cate_draws, "", "")
gh_t <- g_round4(gh_t); names(gh_t)[names(gh_t) == "Estimate"] <- "Diff"
write.csv(gh_t[, c("Covariate","Contrast","Diff","CI_low","CI_high","P_gt_0")],
          file.path(outdir, "gamma_hurdle_hurdle_contrasts.csv"), row.names = FALSE)
g_plot <- gt; names(g_plot)[names(g_plot) == "CATE"] <- "Estimate"
cate_forest_plot(g_plot, "17_ONC_gamma_hurdle_cate_subgroups.png",
                 "Gamma Hurdle: subgroup CATE estimates with 95% CrI")

# ---- 3. Distinct calendar days (both models) -------------------------------
cat("Distinct calendar days, both models...\n")
nonoverlap_file <- file.path(outdir, "New_Data", "analysis",
                             "nonoverlap_model_posterior_draws.rds")
if (file.exists(nonoverlap_file)) {
  nd <- readRDS(nonoverlap_file)
  sub <- rbind(
    build_table(nd$zic$cate, "ZIC-BCF-Smear", "Total non-overlapping severe-AE days"),
    build_table(nd$zic$hurdle_cate, "ZIC-BCF-Smear", "Pr(any non-overlapping severe-AE day)"),
    build_table(nd$gamma$cate, "Gamma hurdle", "Total non-overlapping severe-AE days"),
    build_table(nd$gamma$hurdle_cate, "Gamma hurdle", "Pr(any non-overlapping severe-AE day)"))
  con <- rbind(
    build_contrasts(nd$zic$cate, "ZIC-BCF-Smear", "Total non-overlapping severe-AE days"),
    build_contrasts(nd$zic$hurdle_cate, "ZIC-BCF-Smear", "Pr(any non-overlapping severe-AE day)"),
    build_contrasts(nd$gamma$cate, "Gamma hurdle", "Total non-overlapping severe-AE days"),
    build_contrasts(nd$gamma$hurdle_cate, "Gamma hurdle", "Pr(any non-overlapping severe-AE day)"))
  write.csv(sub, file.path(outdir, "New_Data", "analysis",
                           "nonoverlap_subgroup_model_results.csv"), row.names = FALSE)
  write.csv(con, file.path(outdir, "New_Data", "analysis",
                           "nonoverlap_subgroup_contrasts.csv"), row.names = FALSE)
} else {
  cat("  nonoverlap draws not found; skipping the distinct-day artifacts.\n")
}

cat("Done. Prior SCCHN treatment rows removed from all subgroup artifacts.\n")
