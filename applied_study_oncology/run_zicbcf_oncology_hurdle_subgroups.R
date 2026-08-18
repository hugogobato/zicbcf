# Posterior subgroup analysis for the participation margin of ZIC-BCF-Smear.
# Run run_zicbcf_oncology.R first to create onc_hurdle_draws.rds.

outdir <- "applied_study_oncology"
source(file.path(outdir, "oncology_common.R"))
draw_file <- file.path(outdir, "onc_hurdle_draws.rds")
data_file <- file.path(outdir, "zic_bcf_headneck_analysis_data.csv")
if (!file.exists(draw_file)) stop("Missing hurdle posterior draws: ", draw_file)
if (!file.exists(data_file)) stop("Missing analysis data: ", data_file)

df <- read.csv(data_file, stringsAsFactors = FALSE)
hurdle_cate_draws <- readRDS(draw_file)$hurdle_cate_draws

subgroup_posterior <- function(draws_mat, idx) {
  posterior <- rowMeans(draws_mat[, idx, drop = FALSE])
  c(N = length(idx), CATE = mean(posterior),
    CI_low = unname(quantile(posterior, 0.025)),
    CI_high = unname(quantile(posterior, 0.975)),
    P_gt_0 = mean(posterior > 0))
}

contrast_posterior <- function(draws_mat, idx_a, idx_b) {
  posterior <- rowMeans(draws_mat[, idx_a, drop = FALSE]) -
    rowMeans(draws_mat[, idx_b, drop = FALSE])
  c(Diff = mean(posterior), CI_low = unname(quantile(posterior, 0.025)),
    CI_high = unname(quantile(posterior, 0.975)), P_gt_0 = mean(posterior > 0))
}

group_vars <- onc_subgroup_labels(df)

rows <- list(data.frame(Covariate = "Overall", Level = "All",
                        t(subgroup_posterior(hurdle_cate_draws, seq_len(nrow(df)))),
                        check.names = FALSE))
for (gname in names(group_vars)) {
  lab <- group_vars[[gname]]
  for (level in levels(lab)) {
    idx <- which(lab == level)
    if (length(idx) < ONC_MIN_SUBGROUP_N) next
    rows[[length(rows) + 1]] <- data.frame(Covariate = gname, Level = level,
                                           t(subgroup_posterior(hurdle_cate_draws, idx)),
                                           check.names = FALSE)
  }
}
subgroups <- do.call(rbind, rows)
rownames(subgroups) <- NULL
subgroups[, c("CATE", "CI_low", "CI_high")] <- round(subgroups[, c("CATE", "CI_low", "CI_high")], 4)
subgroups$P_gt_0 <- round(subgroups$P_gt_0, 3)
write.csv(subgroups, file.path(outdir, "onc_hurdle_subgroups.csv"), row.names = FALSE)

rows <- list()
for (gname in names(group_vars)) {
  lab <- group_vars[[gname]]
  eligible <- levels(lab)[table(lab) >= ONC_MIN_SUBGROUP_N]
  if (length(eligible) < 2) next
  pairs <- combn(eligible, 2)
  for (k in seq_len(ncol(pairs))) {
    a <- pairs[1, k]
    b <- pairs[2, k]
    rows[[length(rows) + 1]] <- data.frame(
      Covariate = gname, Contrast = paste(a, "-", b),
      t(contrast_posterior(hurdle_cate_draws, which(lab == a), which(lab == b))),
      check.names = FALSE)
  }
}
contrasts <- do.call(rbind, rows)
rownames(contrasts) <- NULL
contrasts[, c("Diff", "CI_low", "CI_high")] <- round(contrasts[, c("Diff", "CI_low", "CI_high")], 4)
contrasts$P_gt_0 <- round(contrasts$P_gt_0, 3)
write.csv(contrasts, file.path(outdir, "onc_hurdle_contrasts.csv"), row.names = FALSE)
cat("Completed ZIC-BCF-Smear hurdle-margin subgroup analysis.\n")
