# ==========================================================================
# Pool the four final MEPS chain outputs (meps_final_chain{1..4}.rds, written
# by run_meps_final_chain.R, one per Colab runtime or local process) into the
# artifacts the downstream pipeline consumes:
#
#   - meps_zicbcf_posterior_summaries.rds   (input to run_cate_subgroups.R)
#   - meps_unit_level_posterior_means.rds   (marginal figures)
#   - meps_convergence_diagnostics.csv
#   - meps_smearing_{bp_test,bp_draws,residual_scale,sensitivity}.csv
#   - meps_smearing_sensitivity_subgroups.csv (Table: locally smeared subgroups)
#   - meps_smearing_subgroup_relative_change.csv
#
# Usage, from the repository root:
#   Rscript applied_study/combine_meps_final_chains.R
# ==========================================================================

suppressPackageStartupMessages(library(zicbcf))
source("applied_study/meps_common.R")

OUTDIR <- "applied_study"
CHAIN_SEEDS <- seq_len(4L)

chain_files <- file.path(OUTDIR, sprintf("meps_final_chain%d.rds", CHAIN_SEEDS))
if (!all(file.exists(chain_files))) {
  stop("Missing per-chain outputs; run run_meps_final_chain.R 1-4 first.")
}
chains <- lapply(chain_files, readRDS)
dd <- readRDS(file.path(OUTDIR, "meps_final_meta.rds"))
digests <- vapply(chains, function(c) c$domain_digest, character(1))
stopifnot(length(unique(digests)) == 1L, unique(digests) == dd$domain_digest)
stopifnot(length(unique(vapply(chains, function(c) c$n_burn, integer(1)))) == 1L)

df <- dd$df; weights <- dd$weights; domains <- dd$domains
domain_names <- dd$domain_names
total_draws <- sum(vapply(chains, function(c) c$n_draws, integer(1)))
cat("Pooling chains", paste(CHAIN_SEEDS, collapse = ", "),
    "with", total_draws, "total draws (burn-in",
    chains[[1L]]$n_burn, "per chain)\n")

pool <- function(component) do.call(rbind, lapply(chains, function(c) c[[component]]))

# -- monitored functionals -> convergence diagnostics ------------------------
monitored_by_quantity <- lapply(names(chains[[1L]]$monitored), function(q) {
  lapply(chains, function(c) c$monitored[[q]])
})
names(monitored_by_quantity) <- names(chains[[1L]]$monitored)
convergence <- zicbcf_convergence_table(monitored_by_quantity)
write.csv(convergence, file.path(OUTDIR, "meps_convergence_diagnostics.csv"),
          row.names = FALSE)
cat(sprintf("Convergence: max Rhat = %.4f, min ESS = %.0f, max |Geweke z| = %.2f\n",
            max(convergence$rhat, na.rm = TRUE), min(convergence$ess, na.rm = TRUE),
            max(convergence$geweke_z, na.rm = TRUE)))

# -- smearing diagnostics (chain 1) ------------------------------------------
sm1 <- chains[[1L]]$smearing
if (!is.null(sm1)) {
  write.csv(sm1$test, file.path(OUTDIR, "meps_smearing_bp_test.csv"), row.names = FALSE)
  write.csv(sm1$draw_test, file.path(OUTDIR, "meps_smearing_bp_draws.csv"), row.names = FALSE)
  write.csv(sm1$residual_scale, file.path(OUTDIR, "meps_smearing_residual_scale.csv"),
            row.names = FALSE)
  write.csv(sm1$sensitivity, file.path(OUTDIR, "meps_smearing_sensitivity.csv"),
            row.names = FALSE)
  cat(sprintf("Breusch-Pagan p = %.4g; locally smeared ATE differs by %.1f%%\n",
              sm1$test$p_value, 100 * sm1$sensitivity$relative_change))
}

# -- pooled domain draws ------------------------------------------------------
cate_unweighted <- pool("scalar_unweighted")
cate_weighted <- pool("scalar_weighted")
# Per-domain pooled draws: concatenate the four chains' per-draw domain means
# into one 8,000-vector per domain, then bind domains as columns.
hurdle_unweighted <- do.call(cbind, lapply(seq_along(domain_names), function(j)
  unlist(lapply(chains, function(c) c$hurdle_domain[[j]]$unweighted))))
hurdle_weighted <- do.call(cbind, lapply(seq_along(domain_names), function(j)
  unlist(lapply(chains, function(c) c$hurdle_domain[[j]]$weighted))))
stopifnot(dim(hurdle_unweighted) == c(total_draws, length(domain_names)),
          dim(hurdle_weighted) == c(total_draws, length(domain_names)))
colnames(hurdle_unweighted) <- domain_names
colnames(hurdle_weighted) <- domain_names
cate_unit_mean <- colSums(pool("cate_unit_sum")) / total_draws
hurdle_unit_mean <- colSums(pool("hurdle_unit_sum")) / total_draws

# -- design-based variance at the posterior-mean unit effect -----------------
design_variance <- function(unit_mean) {
  vapply(domains, function(d) {
    indicator <- logical(nrow(df)); indicator[d$index] <- TRUE
    meps_design_variance(unit_mean, weights, df$VARSTR, df$VARPSU, indicator)
  }, numeric(1))
}
cate_dv <- design_variance(cate_unit_mean)
hurdle_dv <- design_variance(hurdle_unit_mean)
names(cate_dv) <- domain_names; names(hurdle_dv) <- domain_names

saveRDS(list(
  domains = lapply(domains, function(d) d[c("covariate", "level")]),
  domain_names = domain_names,
  domain_n = vapply(domains, function(d) length(d$index), integer(1)),
  domain_weighted_n = vapply(domains, function(d) sum(weights[d$index]), numeric(1)),
  cate = list(unweighted = cate_unweighted, weighted = cate_weighted),
  hurdle = list(unweighted = hurdle_unweighted, weighted = hurdle_weighted),
  cate_design_variance = cate_dv,
  hurdle_design_variance = hurdle_dv,
  cate_unit_mean = cate_unit_mean,
  hurdle_unit_mean = hurdle_unit_mean,
  n_chains = length(CHAIN_SEEDS), n_burn = chains[[1L]]$n_burn,
  n_sim = chains[[1L]]$n_draws, seeds = CHAIN_SEEDS,
  total_draws = total_draws
), file.path(OUTDIR, "meps_zicbcf_posterior_summaries.rds"))

saveRDS(data.frame(cate_posterior_mean = cate_unit_mean,
                   hurdle_posterior_mean = hurdle_unit_mean,
                   AGE23X = df$AGE23X, FAMINC23 = df$FAMINC23,
                   PERWT23F = df$PERWT23F),
        file.path(OUTDIR, "meps_unit_level_posterior_means.rds"))

# -- locally smeared subgroup table (same construction as before) ------------
local_unwt <- pool("local_unweighted")
local_wt <- pool("local_weighted")
local_unit_mean <- colSums(pool("local_unit_sum")) / total_draws

vapply_domain <- function(unit_values, weighted) {
  vapply(domains, function(d) {
    if (weighted) {
      indicator <- logical(nrow(df)); indicator[d$index] <- TRUE
      meps_design_variance(unit_values, weights, df$VARSTR, df$VARPSU, indicator)
    } else NA_real_
  }, numeric(1))
}
summarize <- function(draws, unit_mean, weighted, smear_label) {
  dv <- if (weighted) vapply_domain(unit_mean, TRUE) else rep(NA_real_, length(domains))
  do.call(rbind, lapply(seq_along(domains), function(j) {
    d <- domains[[j]]
    s <- meps_posterior_summary(draws[, j], dv[j])
    data.frame(Smearing = smear_label,
               Basis = if (weighted) "Population (survey-weighted)" else "Analytic sample",
               Covariate = d$covariate, Level = d$level,
               N = length(d$index), Weighted_N = round(sum(weights[d$index])),
               Estimate = unname(s["Estimate"]),
               CI_low = unname(s["CI_low"]), CI_high = unname(s["CI_high"]),
               P_gt_0 = unname(s["P_gt_0"]),
               Design_SE = unname(s["Design_SE"]),
               Design_CI_low = unname(s["CI_low_design"]),
               Design_CI_high = unname(s["CI_high_design"]))
  }))
}
table17 <- rbind(
  summarize(cate_unweighted, cate_unit_mean, FALSE, "Scalar Duan"),
  summarize(local_unwt, local_unit_mean, FALSE, "Locally smeared"),
  summarize(cate_weighted, cate_unit_mean, TRUE, "Scalar Duan"),
  summarize(local_wt, local_unit_mean, TRUE, "Locally smeared")
)
write.csv(table17, file.path(OUTDIR, "meps_smearing_sensitivity_subgroups.csv"),
          row.names = FALSE)

sens <- do.call(rbind, lapply(seq_along(domain_names), function(j) {
  both <- data.frame(
    Domain = domain_names[j],
    N = length(domains[[j]]$index),
    Scalar_Estimate = mean(cate_unweighted[, j]),
    Local_Estimate = mean(local_unwt[, j]),
    Scalar_CI_low = unname(quantile(cate_unweighted[, j], 0.025)),
    Scalar_CI_high = unname(quantile(cate_unweighted[, j], 0.975)),
    Local_CI_low = unname(quantile(local_unwt[, j], 0.025)),
    Local_CI_high = unname(quantile(local_unwt[, j], 0.975)),
    Relative_change = mean(local_unwt[, j]) / mean(cate_unweighted[, j]) - 1,
    Weighted_Scalar_Estimate = mean(cate_weighted[, j]),
    Weighted_Local_Estimate = mean(local_wt[, j]),
    Weighted_Relative_change = mean(local_wt[, j]) / mean(cate_weighted[, j]) - 1)
  both
}))
write.csv(sens, file.path(OUTDIR, "meps_smearing_subgroup_relative_change.csv"),
          row.names = FALSE)

overall <- which(domain_names == "Overall | All")
cat(sprintf("\nHeadline: sample ATE %.2f [%.2f, %.2f]; population ATE %.2f [%.2f, %.2f]\n",
            mean(cate_unweighted[, overall]),
            quantile(cate_unweighted[, overall], 0.025),
            quantile(cate_unweighted[, overall], 0.975),
            mean(cate_weighted[, overall]),
            quantile(cate_weighted[, overall], 0.025),
            quantile(cate_weighted[, overall], 0.975)))
cat(sprintf("Locally smeared overall: %.2f (%.1f%% vs scalar)\n",
            mean(local_unwt[, overall]),
            100 * (mean(local_unwt[, overall]) / mean(cate_unweighted[, overall]) - 1)))
cat("\nWrote pooled summaries, convergence table, smearing diagnostics and",
    "Table 17 CSVs.\nNext: rerun run_cate_subgroups.R and the figure scripts.\n")
