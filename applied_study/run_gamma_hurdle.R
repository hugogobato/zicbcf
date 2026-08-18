library(zicbcf)

# ==========================================================================
# Gamma-hurdle benchmark for the MEPS dental-insurance study
# (Oganisian et al., 2021).
#
# The cleaning rules, design matrix and subgroup definitions are taken from
# meps_common.R, so the benchmark and ZIC-BCF-Smear are always fit to identical
# data with identical encoding. Any difference between them is therefore
# attributable to the specification and not to the preparation.
#
# As with run_study.R, chains are run and summarized one at a time and only the
# per-draw domain averages are retained, which is sufficient for every reported
# subgroup effect and contrast.
# ==========================================================================

source("applied_study/meps_common.R")

OUTDIR <- "applied_study"
N_CHAINS <- 4L
N_BURN <- 1000L
N_SIM <- 1000L
CHAIN_SEEDS <- seq_len(N_CHAINS)

cat("Loading dataset...\n")
df <- meps_load(file.path(OUTDIR, "h251.csv"))
cat("Cleaned sample size:", nrow(df), "\n")

y <- df$y
z <- df$z
X <- meps_design_matrix(df)
weights <- df$PERWT23F
labels <- meps_subgroup_labels(df)

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

summarize_chain <- function(draws) {
  unweighted <- vapply(domains, function(d) rowMeans(draws[, d$index, drop = FALSE]),
                       numeric(nrow(draws)))
  weighted <- vapply(domains, function(d) meps_weighted_draw_means(draws, weights, d$index),
                     numeric(nrow(draws)))
  colnames(unweighted) <- domain_names
  colnames(weighted) <- domain_names
  list(unweighted = unweighted, weighted = weighted)
}

cate_summaries <- list()
hurdle_summaries <- list()
monitored <- list()
cate_unit_sum <- numeric(nrow(df))
hurdle_unit_sum <- numeric(nrow(df))
total_draws <- 0L

for (chain in seq_along(CHAIN_SEEDS)) {
  cat(sprintf("Fitting Gamma hurdle chain %d of %d...\n", chain, N_CHAINS))
  set.seed(CHAIN_SEEDS[chain])
  fit <- gamma_hurdle(y = y, z = z, x = X, nburn = N_BURN, nsim = N_SIM, nthin = 1)
  hurdle_cate <- fit$p1 - fit$p0

  monitored[[chain]] <- list(
    `ATE (response scale)` = fit$ate,
    `Weighted ATE (response scale)` = meps_weighted_draw_means(fit$cate, weights),
    `Hurdle ATE (probability)` = rowMeans(hurdle_cate)
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

monitored_by_quantity <- lapply(names(monitored[[1L]]), function(q) {
  lapply(monitored, function(m) m[[q]])
})
names(monitored_by_quantity) <- names(monitored[[1L]])
convergence <- zicbcf_convergence_table(monitored_by_quantity)
write.csv(convergence, file.path(OUTDIR, "meps_gamma_convergence_diagnostics.csv"),
          row.names = FALSE)
cat(sprintf("  max Rhat = %.4f, min ESS = %.0f, max |Geweke z| = %.2f\n",
            max(convergence$rhat, na.rm = TRUE), min(convergence$ess, na.rm = TRUE),
            max(convergence$geweke_z, na.rm = TRUE)))

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
), file.path(OUTDIR, "gamma_hurdle_draws.rds"))

overall <- which(domain_names == "Overall | All")
results <- data.frame(
  Estimand = c("Dollar-scale ATE, analytic-sample average",
               "Dollar-scale ATE, survey-weighted population target",
               "Participation-margin ATE, analytic-sample average",
               "Participation-margin ATE, survey-weighted population target"),
  rbind(
    meps_posterior_summary(cate_draws_by_domain$unweighted[, overall]),
    meps_posterior_summary(cate_draws_by_domain$weighted[, overall],
                           cate_design_variance[overall]),
    meps_posterior_summary(hurdle_draws_by_domain$unweighted[, overall]),
    meps_posterior_summary(hurdle_draws_by_domain$weighted[, overall],
                           hurdle_design_variance[overall])
  ),
  row.names = NULL, check.names = FALSE
)
write.csv(results, file.path(OUTDIR, "gamma_hurdle_results.csv"), row.names = FALSE)

png(file.path(OUTDIR, "gamma_hurdle_cate_histogram.png"), width = 800, height = 600)
hist(cate_unit_mean, breaks = 50, col = "skyblue",
     main = "Gamma Hurdle: distribution of unit-level CATE",
     xlab = "CATE (effect of dental insurance on dental expenditure, $)")
dev.off()

png(file.path(OUTDIR, "gamma_hurdle_hurdle_cate_histogram.png"), width = 800, height = 600)
hist(hurdle_unit_mean, breaks = 50, col = "lightgreen",
     main = "Gamma Hurdle: distribution of unit-level participation-margin CATE",
     xlab = "Participation-margin CATE (probability difference)")
dev.off()

cat("\n=== Gamma hurdle average treatment effects ===\n")
print(results, digits = 5, row.names = FALSE)
cat("\nCompleted successfully. Run run_gamma_hurdle_subgroups.R for subgroup tables.\n")
