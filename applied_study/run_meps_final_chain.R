# ==========================================================================
# Final MEPS ZIC-BCF-Smear fit, one chain per process, for Colab or local use.
#
# Review resolution: the Geweke investigation (see
# applied_study_oncology/geweke_investigation/findings.md) showed that the
# flagged family-wise statistics attach to sigma_c because the sigma-tree
# Gibbs coupling equilibrates more slowly than the burn-in previously used.
# Extending burn-in from 1,000 to 5,000 iterations removes every
# family-wise-threatening value without moving posterior summaries beyond
# Monte Carlo noise, so the final models use the longer burn-in.
#
# This script refits ONE chain (argument: 1-4) and saves a single RDS with
# everything downstream inference needs, so the full draw-by-unit matrices
# and the multi-GB raw sampler state never need to leave the runtime:
#
#   - the seven monitored functionals (for the convergence table);
#   - scalar-smear domain draws, unweighted and survey-weighted;
#   - locally smeared domain draws, unweighted and survey-weighted;
#   - unit-level posterior sums of scalar and locally smeared CATEs;
#   - Duan smearing homoskedasticity diagnostics (chain 1 only);
#   - the metadata frame needed to recombine the chains locally.
#
# Usage (from the repository root, after uploading h251.csv to
# applied_study/):
#   Rscript applied_study/run_meps_final_chain.R 1
#   ... one process per chain, then run combine_meps_final_chains.R locally.
# ==========================================================================

suppressPackageStartupMessages(library(zicbcf))
source("applied_study/meps_common.R")

OUTDIR <- "applied_study"
N_BURN <- 5000L
N_SIM <- 2000L
N_CHAINS <- 4L
CHAIN_SEEDS <- seq_len(N_CHAINS)

stage <- commandArgs(trailingOnly = TRUE)[1L]
chain_id <- suppressWarnings(as.integer(stage))
if (is.na(chain_id) || !(chain_id %in% CHAIN_SEEDS)) {
  stop("Usage: run_meps_final_chain.R <chain 1-", N_CHAINS, ">")
}

n_sim_rows <- function(fit) nrow(fit$mu_c)
digest_domains <- function(domains) {
  paste(vapply(domains, function(d) paste(d$covariate, d$level, length(d$index),
                                          sep = "|"), character(1)),
        collapse = ";")
}

# -- data, design, propensity, domains (identical to run_study.R) -----------
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
domain_digest <- digest_domains(domains)
cat(sprintf("[chain %d] data: %d records, %d positive; %d domains\n",
            chain_id, nrow(df), sum(df$y > 0), length(domains)))

if (chain_id == 1L) {
  saveRDS(list(df = df, weights = weights, domains = domains,
               domain_names = domain_names, domain_digest = domain_digest,
               n_burn = N_BURN, n_sim = N_SIM, seeds = CHAIN_SEEDS),
          file.path(OUTDIR, "meps_final_meta.rds"))
}

# -- fit --------------------------------------------------------------------
cat(sprintf("[chain %d] fitting ZIC-BCF-Smear (%d burn-in, %d retained)...\n",
            chain_id, N_BURN, N_SIM))
set.seed(chain_id)
fit <- zicbcf_smear(y = df$y, z = df$z, x_control = X, x_moderate = X,
                    pihat = pihat, nburn = N_BURN, nsim = N_SIM)

# -- homoskedasticity diagnostics (chain 1 only, as in run_study.R) ---------
smearing <- NULL
if (chain_id == 1L) {
  cat("[chain 1] testing the homoskedasticity assumption behind Duan's smearing...\n")
  smearing <- zicbcf_smearing_diagnostics(fit, y = df$y, z = df$z, x = X)
}

# -- monitored functionals (convergence table) ------------------------------
hurdle_cate <- pnorm(fit$mu_b + fit$tau_b) - pnorm(fit$mu_b)
monitored <- list(
  `ATE (response scale)` = fit$ate,
  `Weighted ATE (response scale)` = meps_weighted_draw_means(fit$cate, weights),
  `Hurdle ATE (probability)` = rowMeans(hurdle_cate),
  `Duan smearing factor` = fit$smearing_factors,
  `sigma_c (log-scale SD)` = as.numeric(fit$sigma_c),
  `Mean mu_c (log-scale prognostic)` = rowMeans(fit$mu_c),
  `Mean tau_c (log-scale treatment)` = rowMeans(fit$tau_c)
)

# -- locally smeared response-scale CATEs, draw by draw ---------------------
active <- which(df$y > 0)
log_y_active <- log(df$y[active])
z_active <- as.numeric(df$z[active])
n_units <- nrow(X)
n_bins <- 5L

mu_b <- fit$mu_b; tau_b <- fit$tau_b
mu_c <- fit$mu_c; tau_c <- fit$tau_c
cate_scalar <- fit$cate
smearing_factors <- as.numeric(fit$smearing_factors)
ate_chain <- fit$ate
rm(fit); invisible(gc(verbose = FALSE))

fitted_active <- mu_c[, active, drop = FALSE] +
  matrix(z_active, nrow = N_SIM, ncol = length(active), byrow = TRUE) *
  tau_c[, active, drop = FALSE]
fitted_mean <- colMeans(fitted_active)
breaks <- stats::quantile(fitted_mean, probs = seq(0, 1, length.out = n_bins + 1L))
bin_edges <- unique(breaks)
n_strata <- length(bin_edges) - 1L
assign_bin <- function(v) {
  clamped <- pmin(pmax(v, bin_edges[1L]), bin_edges[length(bin_edges)])
  as.integer(cut(clamped, breaks = bin_edges, include.lowest = TRUE,
                 labels = FALSE))
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
local_unwt <- matrix(NA_real_, N_SIM, length(domains))
local_wt <- matrix(NA_real_, N_SIM, length(domains))
local_unit_sum <- numeric(n_units)

t0 <- Sys.time()
for (s in seq_len(N_SIM)) {
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
  local_unwt[s, ] <- vapply(domain_indices, function(idx) mean(cate_local_s[idx]),
                            numeric(1))
  local_wt[s, ] <- vapply(domain_indices, function(idx) {
    sw <- sum(weights[idx]); if (sw <= 0) return(NA_real_)
    sum(weights[idx] * cate_local_s[idx]) / sw
  }, numeric(1))
  if (s %% 500L == 0L) {
    cat(sprintf("[chain %d] locally smeared draw %d of %d (%.1f min elapsed)\n",
                chain_id, s, N_SIM,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    invisible(gc(verbose = FALSE))
  }
}
rm(cate_local_s); invisible(gc(verbose = FALSE))

# -- scalar domain summaries -------------------------------------------------
scalar_unwt <- vapply(domain_indices, function(idx)
  rowMeans(cate_scalar[, idx, drop = FALSE]), numeric(N_SIM))
scalar_wt <- vapply(domain_indices, function(idx)
  meps_weighted_draw_means(cate_scalar, weights, idx), numeric(N_SIM))

hurdle_summaries <- lapply(domain_indices, function(idx)
  list(unweighted = rowMeans(hurdle_cate[, idx, drop = FALSE]),
       weighted = meps_weighted_draw_means(hurdle_cate, weights, idx)))

saveRDS(list(
  chain_id = chain_id,
  monitored = monitored,
  scalar_unweighted = scalar_unwt, scalar_weighted = scalar_wt,
  local_unweighted = local_unwt, local_weighted = local_wt,
  cate_unit_sum = colSums(cate_scalar), local_unit_sum = local_unit_sum,
  hurdle_domain = hurdle_summaries,
  hurdle_unit_sum = colSums(hurdle_cate),
  smearing = smearing,
  smearing_factor_scalar = smearing_factors,
  residual_scale = residual_scale,
  bin_edges = bin_edges,
  ate_chain = ate_chain,
  n_draws = N_SIM,
  n_burn = N_BURN,
  domain_digest = domain_digest
), file.path(OUTDIR, sprintf("meps_final_chain%d.rds", chain_id)))

overall <- which(domain_names == "Overall | All")
cat(sprintf("[chain %d] done: scalar ATE %.2f, locally smeared ATE %.2f (%.1f%%)\n",
            chain_id, mean(ate_chain), mean(local_unwt[, overall]),
            100 * (mean(local_unwt[, overall]) / mean(ate_chain) - 1)))
