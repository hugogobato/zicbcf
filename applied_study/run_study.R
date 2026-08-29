library(zicbcf)

# ==========================================================================
# MEPS 2023 applied study: the causal effect of dental-insurance coverage on
# annual dental expenditure. This script performs the ZIC-BCF-Smear fit and
# writes the posterior summaries that run_cate_subgroups.R turns into tables
# and figures.
#
# The cleaning rules, design matrix, subgroup definitions and survey-design
# handling live in meps_common.R; read the header of that file for what changed
# relative to the first version of this analysis.
#
# WHAT IS SAVED, AND WHY IT IS NOT THE FULL DRAW MATRIX. Every reported
# quantity (subgroup effects, between-subgroup contrasts, and their weighted
# counterparts) is an average of unit-level effects taken WITHIN a posterior
# draw. Those per-draw averages are therefore sufficient for all downstream
# inference, and they are four orders of magnitude smaller than the
# draws-by-units matrix. Unit-level posterior means are saved separately for
# the marginal figures. Chains are run one at a time and discarded after
# summarization, which keeps peak memory bounded on an 18,763-record sample.
# ==========================================================================

source("applied_study/meps_common.R")

OUTDIR <- "applied_study"
N_CHAINS <- 4L
N_BURN <- 5000L  # raised from 1,000 after the Geweke investigation; final models
N_SIM <- 2000L
CHAIN_SEEDS <- seq_len(N_CHAINS)

cat("Loading dataset...\n")
df <- meps_load(file.path(OUTDIR, "h251.csv"))
cat("Cleaned sample size:", nrow(df), "\n")

y <- df$y
z <- df$z
X <- meps_design_matrix(df)
weights <- df$PERWT23F
labels <- meps_subgroup_labels(df)

cat("Design matrix columns:", paste(colnames(X), collapse = ", "), "\n")
cat(sprintf("Zeros in outcome: %d (%.1f%%); positive median $%.0f, max $%.0f\n",
            sum(y == 0), 100 * mean(y == 0), median(y[y > 0]), max(y)))
cat(sprintf("Survey design: %d strata, %d primary sampling units, %d records with zero weight\n",
            length(unique(df$VARSTR)), nrow(unique(df[, c("VARSTR", "VARPSU")])),
            sum(weights == 0)))
cat(sprintf("Weighted population represented: %s persons\n",
            format(round(sum(weights)), big.mark = ",")))

cat("Estimating propensity scores...\n")
ps_model <- glm(z ~ ., data = as.data.frame(X), family = binomial())
pihat <- predict(ps_model, type = "response")
cat("pihat range: [", round(min(pihat), 3), ",", round(max(pihat), 3), "]\n")

# ---------------------------------------------------------------------------
# Fixed list of reported domains: the overall sample plus every subgroup level
# that meets the minimum size.
# ---------------------------------------------------------------------------
domains <- list(list(covariate = "Overall", level = "All",
                     index = seq_len(nrow(df))))
for (group_name in names(labels)) {
  label <- labels[[group_name]]
  for (level in levels(label)) {
    index <- which(label == level)
    if (length(index) < MEPS_MIN_SUBGROUP_N) next
    domains[[length(domains) + 1L]] <- list(covariate = group_name,
                                            level = level, index = index)
  }
}
domain_names <- vapply(domains, function(d) paste(d$covariate, d$level, sep = " | "),
                       character(1))
cat("Reported domains:", length(domains), "\n")

summarize_chain <- function(cate_draws) {
  unweighted <- vapply(domains, function(d) {
    rowMeans(cate_draws[, d$index, drop = FALSE])
  }, numeric(nrow(cate_draws)))
  weighted <- vapply(domains, function(d) {
    meps_weighted_draw_means(cate_draws, weights, d$index)
  }, numeric(nrow(cate_draws)))
  colnames(unweighted) <- domain_names
  colnames(weighted) <- domain_names
  list(unweighted = unweighted, weighted = weighted)
}

# ---------------------------------------------------------------------------
# Chains, run and summarized one at a time.
# ---------------------------------------------------------------------------
cate_summaries <- list()
hurdle_summaries <- list()
monitored <- list()
cate_unit_sum <- numeric(nrow(df))
hurdle_unit_sum <- numeric(nrow(df))
total_draws <- 0L
smearing <- NULL

for (chain in seq_along(CHAIN_SEEDS)) {
  seed <- CHAIN_SEEDS[chain]
  cat(sprintf("Fitting ZIC-BCF-Smear chain %d of %d...\n", chain, N_CHAINS))
  set.seed(seed)
  fit <- zicbcf_smear(y = y, z = z, x_control = X, x_moderate = X,
                      pihat = pihat, nburn = N_BURN, nsim = N_SIM)

  if (chain == 1L) {
    cat("  testing the homoskedasticity assumption behind Duan's smearing...\n")
    smearing <- zicbcf_smearing_diagnostics(fit, y = y, z = z, x = X)
  }

  hurdle_cate <- pnorm(fit$mu_b + fit$tau_b) - pnorm(fit$mu_b)
  monitored[[chain]] <- list(
    `ATE (response scale)` = fit$ate,
    `Weighted ATE (response scale)` = meps_weighted_draw_means(fit$cate, weights),
    `Hurdle ATE (probability)` = rowMeans(hurdle_cate),
    `Duan smearing factor` = fit$smearing_factors,
    `sigma_c (log-scale SD)` = as.numeric(fit$sigma_c),
    `Mean mu_c (log-scale prognostic)` = rowMeans(fit$mu_c),
    `Mean tau_c (log-scale treatment)` = rowMeans(fit$tau_c)
  )

  cate_summaries[[chain]] <- summarize_chain(fit$cate)
  hurdle_summaries[[chain]] <- summarize_chain(hurdle_cate)
  cate_unit_sum <- cate_unit_sum + colSums(fit$cate)
  hurdle_unit_sum <- hurdle_unit_sum + colSums(hurdle_cate)
  total_draws <- total_draws + nrow(fit$cate)

  rm(fit, hurdle_cate)
  invisible(gc(verbose = FALSE))
}

pool <- function(summaries, component) {
  do.call(rbind, lapply(summaries, function(s) s[[component]]))
}
cate_draws_by_domain <- list(unweighted = pool(cate_summaries, "unweighted"),
                             weighted = pool(cate_summaries, "weighted"))
hurdle_draws_by_domain <- list(unweighted = pool(hurdle_summaries, "unweighted"),
                               weighted = pool(hurdle_summaries, "weighted"))
cate_unit_mean <- cate_unit_sum / total_draws
hurdle_unit_mean <- hurdle_unit_sum / total_draws

# ---------------------------------------------------------------------------
# Convergence diagnostics
# ---------------------------------------------------------------------------
cat("Computing convergence diagnostics...\n")
monitored_by_quantity <- lapply(names(monitored[[1L]]), function(q) {
  lapply(monitored, function(m) m[[q]])
})
names(monitored_by_quantity) <- names(monitored[[1L]])
convergence <- zicbcf_convergence_table(monitored_by_quantity)
write.csv(convergence, file.path(OUTDIR, "meps_convergence_diagnostics.csv"),
          row.names = FALSE)
cat(sprintf("  max Rhat = %.4f, min ESS = %.0f, max |Geweke z| = %.2f\n",
            max(convergence$rhat, na.rm = TRUE), min(convergence$ess, na.rm = TRUE),
            max(convergence$geweke_z, na.rm = TRUE)))

write.csv(smearing$test, file.path(OUTDIR, "meps_smearing_bp_test.csv"), row.names = FALSE)
write.csv(smearing$draw_test, file.path(OUTDIR, "meps_smearing_bp_draws.csv"), row.names = FALSE)
write.csv(smearing$residual_scale, file.path(OUTDIR, "meps_smearing_residual_scale.csv"),
          row.names = FALSE)
write.csv(smearing$sensitivity, file.path(OUTDIR, "meps_smearing_sensitivity.csv"),
          row.names = FALSE)
cat(sprintf("  Breusch-Pagan p = %.4g; locally smeared ATE differs by %.1f%%\n",
            smearing$test$p_value, 100 * smearing$sensitivity$relative_change))

# ---------------------------------------------------------------------------
# Design-based variance of every reported domain mean, evaluated at the
# posterior-mean unit-level effect.
# ---------------------------------------------------------------------------
cat("Computing design-based variances...\n")
design_variance <- function(unit_mean) {
  vapply(domains, function(d) {
    domain_indicator <- logical(nrow(df))
    domain_indicator[d$index] <- TRUE
    meps_design_variance(unit_mean, weights, df$VARSTR, df$VARPSU, domain_indicator)
  }, numeric(1))
}
cate_design_variance <- design_variance(cate_unit_mean)
hurdle_design_variance <- design_variance(hurdle_unit_mean)
names(cate_design_variance) <- domain_names
names(hurdle_design_variance) <- domain_names

# ---------------------------------------------------------------------------
# Save everything downstream scripts need
# ---------------------------------------------------------------------------
saveRDS(list(
  domains = lapply(domains, function(d) d[c("covariate", "level")]),
  domain_names = domain_names,
  domain_n = vapply(domains, function(d) length(d$index), integer(1)),
  domain_weighted_n = vapply(domains, function(d) sum(weights[d$index]), numeric(1)),
  cate = cate_draws_by_domain,
  hurdle = hurdle_draws_by_domain,
  cate_design_variance = cate_design_variance,
  hurdle_design_variance = hurdle_design_variance,
  cate_unit_mean = cate_unit_mean,
  hurdle_unit_mean = hurdle_unit_mean,
  n_chains = N_CHAINS, n_burn = N_BURN, n_sim = N_SIM, seeds = CHAIN_SEEDS,
  total_draws = total_draws
), file.path(OUTDIR, "meps_zicbcf_posterior_summaries.rds"))

saveRDS(data.frame(cate_posterior_mean = cate_unit_mean,
                   hurdle_posterior_mean = hurdle_unit_mean,
                   AGE23X = df$AGE23X, FAMINC23 = df$FAMINC23,
                   PERWT23F = df$PERWT23F),
        file.path(OUTDIR, "meps_unit_level_posterior_means.rds"))

# ---------------------------------------------------------------------------
# Overall average treatment effect, sample average and population target
# ---------------------------------------------------------------------------
overall <- which(domain_names == "Overall | All")
sample_average <- meps_posterior_summary(cate_draws_by_domain$unweighted[, overall])
population <- meps_posterior_summary(cate_draws_by_domain$weighted[, overall],
                                     cate_design_variance[overall])
hurdle_sample <- meps_posterior_summary(hurdle_draws_by_domain$unweighted[, overall])
hurdle_population <- meps_posterior_summary(hurdle_draws_by_domain$weighted[, overall],
                                            hurdle_design_variance[overall])

ate_results <- data.frame(
  Estimand = c("Dollar-scale ATE, analytic-sample average",
               "Dollar-scale ATE, survey-weighted population target",
               "Participation-margin ATE, analytic-sample average",
               "Participation-margin ATE, survey-weighted population target"),
  rbind(sample_average, population, hurdle_sample, hurdle_population),
  row.names = NULL, check.names = FALSE
)
write.csv(ate_results, file.path(OUTDIR, "ate_results.csv"), row.names = FALSE)

png(file.path(OUTDIR, "cate_histogram.png"), width = 800, height = 600)
hist(cate_unit_mean, breaks = 50, col = "skyblue",
     main = "Distribution of Conditional Average Treatment Effects (CATE)",
     xlab = "CATE (impact of dental insurance on dental expenditure, $)")
dev.off()

cat("\n=== Average treatment effects ===\n")
print(ate_results, digits = 5, row.names = FALSE)
cat("\nCompleted successfully. Run run_cate_subgroups.R for the subgroup tables.\n")
