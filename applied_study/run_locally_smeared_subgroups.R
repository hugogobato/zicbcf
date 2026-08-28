# ==========================================================================
# Locally smeared sensitivity analysis for the MEPS subgroup table.
#
# Review point: the locally smeared sensitivity analysis existed only for the
# population-average effect, but the smearing factor rises monotonically
# across strata of the fitted log-scale mean, so subgroup ordering (what the
# disparity gradient rests on) is exposed to residual-scale misscaling in a
# way that an average over all units partially cancels. This script
# recomputes every reported subgroup effect with the single scalar smearing
# factor replaced, within each posterior draw, by stratum-specific factors
# over quantiles of the fitted log-scale mean, exactly as
# zicbcf_smearing_diagnostics() does for the average.
#
# The full per-draw log-scale surfaces were discarded after summarization by
# run_study.R, so this script refits the same four chains with the same seeds
# (1 through 4), burn-in (1,000) and retained draws (2,000) and saves what
# downstream inference needs:
#
#   - scalar-smear domain draws, unweighted and survey-weighted;
#   - locally smeared domain draws, unweighted and survey-weighted;
#   - unit-level posterior means of both smear types, for design variance.
#
# Usage: run one chain per process so the four fits use separate cores, then
# pool them.
#   Rscript applied_study/run_locally_smeared_subgroups.R 1
#   Rscript applied_study/run_locally_smeared_subgroups.R 2
#   Rscript applied_study/run_locally_smeared_subgroups.R 3
#   Rscript applied_study/run_locally_smeared_subgroups.R 4
#   Rscript applied_study/run_locally_smeared_subgroups.R combine
#
# Run from the repository root. Requires the zicbcf package to be installed
# (R CMD INSTALL . at the repository root) and applied_study/h251.csv to be
# present (fetch with download_meps_h251.R).
# ==========================================================================

library(zicbcf)
source("applied_study/meps_common.R")

OUTDIR <- "applied_study"
N_CHAINS <- 4L
N_BURN <- 1000L
N_SIM <- 2000L
CHAIN_SEEDS <- seq_len(N_CHAINS)

n_sim_rows <- function(fit) nrow(fit$mu_c)
digest_domains <- function(domains) {
  paste(vapply(domains, function(d) paste(d$covariate, d$level, length(d$index),
                                          sep = "|"), character(1)),
        collapse = ";")
}

stage <- commandArgs(trailingOnly = TRUE)[1L]
if (is.na(stage)) stop("Usage: run_locally_smeared_subgroups.R <chain 1-4 | combine>")

domains_and_data <- function() {
  df <- meps_load(file.path(OUTDIR, "h251.csv"))
  X <- meps_design_matrix(df)
  ps_model <- glm(z ~ ., data = cbind(z = df$z, as.data.frame(X)), family = binomial())
  pihat <- predict(ps_model, type = "response")
  weights <- df$PERWT23F
  labels <- meps_subgroup_labels(df)
  domains <- list(list(covariate = "Overall", level = "All",
                       index = seq_len(nrow(df))))
  for (group_name in names(labels)) {
    label <- labels[[group_name]]
    for (level in levels(label)) {
      index <- which(label == level)
      if (length(index) < MEPS_MIN_SUBGROUP_N) next
      domains[[length(domains) + 1L]] <-
        list(covariate = group_name, level = level, index = index)
    }
  }
  domain_names <- vapply(domains, function(d) paste(d$covariate, d$level, sep = " | "),
                         character(1))
  list(df = df, X = X, pihat = pihat, weights = weights,
       domains = domains, domain_names = domain_names)
}

if (stage == "combine") {
  chain_files <- file.path(OUTDIR, sprintf("smearing_subgroup_chain%d.rds", CHAIN_SEEDS))
  if (!all(file.exists(chain_files))) {
    stop("Missing per-chain outputs; run chains 1-", N_CHAINS, " first.")
  }
  cat("Pools chains:", paste(CHAIN_SEEDS, collapse = ", "), "\n")
  chains <- lapply(chain_files, readRDS)
  dd <- readRDS(file.path(OUTDIR, "smearing_subgroup_meta.rds"))
  digests <- vapply(chains, function(c) c$domain_digest, character(1))
  stopifnot(length(unique(digests)) == 1L, unique(digests) == dd$domain_digest)

  pool <- function(component) do.call(rbind, lapply(chains, function(c) c[[component]]))
  scalar_unwt <- pool("scalar_unweighted")
  scalar_wt <- pool("scalar_weighted")
  local_unwt <- pool("local_unweighted")
  local_wt <- pool("local_weighted")
  local_unit_mean <- colSums(pool("local_unit_sum")) / sum(vapply(chains, function(c) c$n_draws, numeric(1)))
  scalar_unit_mean <- colSums(pool("scalar_unit_sum")) / sum(vapply(chains, function(c) c$n_draws, numeric(1)))
  domains <- dd$domains; domain_names <- dd$domain_names
  weights <- dd$weights; df <- dd$df

  vapply_domain <- function(unit_values, weighted) {
    vapply(domains, function(d) {
      if (weighted) {
        meps_design_variance(unit_values, weights, df$VARSTR, df$VARPSU,
                             seq_along(unit_values) %in% d$index)
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
  out <- rbind(
    summarize(scalar_unwt, scalar_unit_mean, FALSE, "Scalar Duan"),
    summarize(local_unwt, local_unit_mean, FALSE, "Locally smeared"),
    summarize(scalar_wt, scalar_unit_mean, TRUE, "Scalar Duan"),
    summarize(local_wt, local_unit_mean, TRUE, "Locally smeared")
  )
  write.csv(out, file.path(OUTDIR, "meps_smearing_sensitivity_subgroups.csv"),
            row.names = FALSE)

  # Same-domain scalar versus local comparison, one row per domain and basis.
  sens <- do.call(rbind, lapply(seq_along(domain_names), function(j) {
    both <- data.frame(
      Domain = domain_names[j],
      N = length(domains[[j]]$index),
      Scalar_Estimate = mean(scalar_unwt[, j]), Local_Estimate = mean(local_unwt[, j]),
      Scalar_CI_low = unname(quantile(scalar_unwt[, j], 0.025)),
      Scalar_CI_high = unname(quantile(scalar_unwt[, j], 0.975)),
      Local_CI_low = unname(quantile(local_unwt[, j], 0.025)),
      Local_CI_high = unname(quantile(local_unwt[, j], 0.975)),
      Relative_change = mean(local_unwt[, j]) / mean(scalar_unwt[, j]) - 1)
    if (!is.null(scalar_wt)) {
      wl <- c(Weighted_Scalar_Estimate = mean(scalar_wt[, j]),
              Weighted_Local_Estimate = mean(local_wt[, j]),
              Weighted_Relative_change =
                mean(local_wt[, j]) / mean(scalar_wt[, j]) - 1)
      both <- cbind(both, t(wl))
    }
    both
  }))
  write.csv(sens, file.path(OUTDIR, "meps_smearing_subgroup_relative_change.csv"),
            row.names = FALSE)
  cat("\n=== Scalar vs locally smeared subgroup estimates (analytic sample) ===\n")
  print(transform(sens, Relative_change = round(Relative_change, 4),
                  Weighted_Relative_change = round(Weighted_Relative_change, 4)),
        digits = 5, row.names = FALSE)
  cat("\nWrote meps_smearing_sensitivity_subgroups.csv and",
      "meps_smearing_subgroup_relative_change.csv\n")
} else {
  chain_id <- as.integer(stage)
  if (!(chain_id %in% CHAIN_SEEDS)) stop("chain must be one of ", paste(CHAIN_SEEDS, collapse = ", "))
  set.seed(chain_id)

  env <- domains_and_data()
  df <- env$df; X <- env$X; pihat <- env$pihat
  weights <- env$weights; domains <- env$domains; domain_names <- env$domain_names

  if (chain_id == 1L) {
    saveRDS(list(df = df, weights = weights, domains = domains,
                 domain_names = domain_names,
                 domain_digest = digest_domains(domains)),
            file.path(OUTDIR, "smearing_subgroup_meta.rds"))
  }

  cat(sprintf("[chain %d] fitting ZIC-BCF-Smear (%d burn-in, %d retained)...\n",
              chain_id, N_BURN, N_SIM))
  fit <- zicbcf_smear(y = df$y, z = df$z, x_control = X, x_moderate = X,
                      pihat = pihat, nburn = N_BURN, nsim = N_SIM)

  # The returned object carries several draw-by-unit matrices plus the raw
  # sampler state (multi-GB). Extract what is needed and drop the rest before
  # the per-draw loop so peak memory stays bounded.
  n_sim <- N_SIM
  mu_b <- fit$mu_b; tau_b <- fit$tau_b
  mu_c <- fit$mu_c; tau_c <- fit$tau_c
  cate_scalar <- fit$cate
  smearing_factors <- as.numeric(fit$smearing_factors)
  ate_chain <- fit$ate
  rm(fit); invisible(gc(verbose = FALSE))
  stopifnot(nrow(mu_c) == n_sim, nrow(cate_scalar) == n_sim)

  # --- locally smeared response-scale CATEs, draw by draw -------------------
  active <- which(df$y > 0)
  log_y_active <- log(df$y[active])
  z_active <- as.numeric(df$z[active])
  n_units <- nrow(X)
  n_bins <- 5L

  fitted_active <- mu_c[, active, drop = FALSE] +
    matrix(z_active, nrow = n_sim, ncol = length(active), byrow = TRUE) *
    tau_c[, active, drop = FALSE]
  fitted_mean <- colMeans(fitted_active)
  breaks <- stats::quantile(fitted_mean, probs = seq(0, 1, length.out = n_bins + 1L))
  bin_edges <- unique(breaks)
  n_strata <- length(bin_edges) - 1L
  assign_bin <- function(v) {
    clamped <- pmin(pmax(v, bin_edges[1L]), bin_edges[length(bin_edges)])
    as.integer(cut(clamped, breaks = bin_edges, include.lowest = TRUE, labels = FALSE))
  }
  resid_mean_active <- log_y_active - fitted_mean
  residual_scale <- do.call(rbind, lapply(sort(unique(assign_bin(fitted_mean))), function(b) {
    keep <- assign_bin(fitted_mean) == b
    data.frame(stratum = b, n = sum(keep),
               fitted_mean_midpoint = stats::median(fitted_mean[keep]),
               residual_sd = stats::sd(resid_mean_active[keep]),
               smearing_factor = mean(exp(resid_mean_active[keep])))
  }))

  domain_indices <- lapply(domains, `[[`, "index")
  local_unwt <- matrix(NA_real_, n_sim, length(domains))
  local_wt <- matrix(NA_real_, n_sim, length(domains))
  local_unit_sum <- numeric(n_units)

  t0 <- Sys.time()
  for (s in seq_len(n_sim)) {
    resid_s <- log_y_active - fitted_active[s, ]
    bin_s <- assign_bin(fitted_active[s, ])
    phi_stratum <- rep(smearing_factors[s], n_strata)
    observed <- tapply(exp(resid_s), factor(bin_s, levels = seq_len(n_strata)), mean)
    phi_stratum[!is.na(observed)] <- observed[!is.na(observed)]

    mu_c_s <- mu_c[s, ]; tau_c_s <- tau_c[s, ]
    phi0 <- phi_stratum[assign_bin(mu_c_s)]
    phi1 <- phi_stratum[assign_bin(mu_c_s + tau_c_s)]
    p0 <- pnorm(mu_b[s, ]); p1 <- pnorm(mu_b[s, ] + tau_b[s, ])
    cate_local_s <- p1 * exp(mu_c_s + tau_c_s) * phi1 -
      p0 * exp(mu_c_s) * phi0
    local_unit_sum <- local_unit_sum + cate_local_s
    local_unwt[s, ] <- vapply(domain_indices, function(idx) mean(cate_local_s[idx]), numeric(1))
    w <- weights
    local_wt[s, ] <- vapply(domain_indices, function(idx) {
      sw <- sum(w[idx]); if (sw <= 0) return(NA_real_)
      sum(w[idx] * cate_local_s[idx]) / sw
    }, numeric(1))
    if (s %% 500L == 0L) {
      cat(sprintf("[chain %d] locally smeared draw %d of %d (%.1f min elapsed)\n",
                  chain_id, s, n_sim, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      invisible(gc(verbose = FALSE))
    }
  }
  rm(cate_local_s); invisible(gc(verbose = FALSE))

  # draw x domain matrices: row s is posterior draw s, column j is domain j.
  scalar_unwt <- vapply(domain_indices, function(idx)
    rowMeans(cate_scalar[, idx, drop = FALSE]), numeric(n_sim))
  scalar_wt <- vapply(domain_indices, function(idx)
    meps_weighted_draw_means(cate_scalar, weights, idx), numeric(n_sim))

  saveRDS(list(
    scalar_unweighted = scalar_unwt, scalar_weighted = scalar_wt,
    local_unweighted = local_unwt, local_weighted = local_wt,
    scalar_unit_sum = colSums(cate_scalar), local_unit_sum = local_unit_sum,
    smearing_factor_scalar = smearing_factors,
    residual_scale = residual_scale,
    bin_edges = bin_edges,
    ate_chain = ate_chain,
    n_draws = n_sim,
    domain_digest = digest_domains(domains)
  ), file.path(OUTDIR, sprintf("smearing_subgroup_chain%d.rds", chain_id)))

  overall <- which(domain_names == "Overall | All")
  cat(sprintf("[chain %d] done: scalar ATE %.2f, locally smeared ATE %.2f (%.1f%%)\n",
              chain_id, mean(ate_chain), mean(local_unwt[, overall]),
              100 * (mean(local_unwt[, overall]) / mean(ate_chain) - 1)))
}
