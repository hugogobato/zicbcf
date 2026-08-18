################################################################################
##  Convergence and smearing diagnostics for the simulation grid.
##
##  The simulation study reports posterior summaries from single chains. This
##  script supplies the diagnostics behind them: for each of the six
##  data-generating processes it runs four independent chains on each of several
##  replicates and reports the Gelman-Rubin statistic, the effective sample
##  size, and the Geweke statistic for the average treatment effect, the
##  log-scale residual standard deviation, Duan's smearing factor, and ten
##  unit-level conditional average treatment effects spanning the estimated
##  effect surface.
##
##  It also runs the studentized Breusch-Pagan test of the homoskedasticity
##  assumption that the scalar smearing factor rests on, together with the
##  locally smeared sensitivity analysis, so that the consequences of a
##  violation are measured rather than assumed. DGP A and its nonlinear
##  counterpart have homoskedastic log-scale errors by construction and
##  therefore act as negative controls; DGP B (Gamma) and DGP C (Tweedie) do
##  not, and act as positive ones.
##
##  Usage:  Rscript simulation_studies/run_convergence_diagnostics.R [n_seeds]
################################################################################

suppressPackageStartupMessages({
  library(zicbcf)
  library(dplyr)
})

RESULTS_DIR <- "simulation_studies/results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
N_SEEDS <- if (length(args) >= 1L) as.integer(args[1L]) else 5L
N_CHAINS <- 4L
NBURN <- 1000L
NSIM <- 1000L
N <- 500L
P <- 5L

source("simulation_studies/dgps_linear.R")
source("simulation_studies/nonlinear_simulation_studies/nonlinear_dgps.R")

generators <- list(
  `Linear DGP A (log-normal hurdle)` = function(seed) generate_dgp_a(seed),
  `Linear DGP B (Gamma hurdle)` = function(seed) generate_dgp_b(seed),
  `Linear DGP C (Tweedie)` = function(seed) generate_dgp_c(seed),
  `Nonlinear DGP A (log-normal hurdle)` = function(seed) generate_dgp_a_nl(N, P, seed),
  `Nonlinear DGP B (Gamma hurdle)` = function(seed) generate_dgp_b_nl(N, P, seed),
  `Nonlinear DGP C (Tweedie)` = function(seed) generate_dgp_c_nl(N, P, seed)
)

convergence_rows <- list()
smearing_rows <- list()

for (dgp_name in names(generators)) {
  for (seed in seq_len(N_SEEDS)) {
    data <- generators[[dgp_name]](seed)
    cat(sprintf("%s, replicate %d: fitting %d chains...\n", dgp_name, seed, N_CHAINS))

    chains <- lapply(seq_len(N_CHAINS), function(chain) {
      set.seed(1000L * seed + chain)
      zicbcf_smear(y = data$y, z = data$z, x_control = data$x,
                   x_moderate = data$x, pihat = data$pihat,
                   nburn = NBURN, nsim = NSIM)
    })

    diagnostics <- zicbcf_convergence(chains, n_cate_units = 10L)
    diagnostics$dgp <- dgp_name
    diagnostics$replicate <- seed
    convergence_rows[[length(convergence_rows) + 1L]] <- diagnostics

    smearing <- zicbcf_smearing_diagnostics(chains[[1L]], y = data$y, z = data$z,
                                            x = data$x)
    smearing_rows[[length(smearing_rows) + 1L]] <- data.frame(
      dgp = dgp_name,
      replicate = seed,
      bp_statistic = smearing$test$statistic,
      bp_df = smearing$test$df,
      bp_p_value = smearing$test$p_value,
      draw_median_p = smearing$draw_test$median_p,
      draw_fraction_p_below_0.05 = smearing$draw_test$p_below_0.05,
      residual_sd_ratio = max(smearing$residual_scale$residual_sd) /
        min(smearing$residual_scale$residual_sd),
      scalar_smearing_ate = smearing$sensitivity$scalar_smearing_ate,
      local_smearing_ate = smearing$sensitivity$local_smearing_ate,
      relative_change = smearing$sensitivity$relative_change,
      true_ate = data$true_ate,
      stringsAsFactors = FALSE
    )

    rm(chains)
    invisible(gc(verbose = FALSE))
  }
}

convergence <- bind_rows(convergence_rows)
smearing <- bind_rows(smearing_rows)
write.csv(convergence, file.path(RESULTS_DIR, "simulation_convergence_diagnostics.csv"),
          row.names = FALSE)
write.csv(smearing, file.path(RESULTS_DIR, "simulation_smearing_diagnostics.csv"),
          row.names = FALSE)

convergence_summary <- convergence %>%
  mutate(kind = if_else(grepl("^CATE\\[", quantity), "Unit-level CATE", "Global summary")) %>%
  group_by(dgp, kind) %>%
  summarise(
    monitored = n(),
    median_rhat = median(rhat, na.rm = TRUE),
    max_rhat = max(rhat, na.rm = TRUE),
    p95_rhat = quantile(rhat, 0.95, na.rm = TRUE),
    median_ess = median(ess, na.rm = TRUE),
    min_ess = min(ess, na.rm = TRUE),
    max_abs_geweke = max(geweke_z, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(convergence_summary,
          file.path(RESULTS_DIR, "simulation_convergence_summary.csv"), row.names = FALSE)

smearing_summary <- smearing %>%
  group_by(dgp) %>%
  summarise(
    replicates = n(),
    median_bp_p = median(bp_p_value),
    rejected_at_0.05 = mean(bp_p_value < 0.05),
    median_residual_sd_ratio = median(residual_sd_ratio),
    median_relative_change = median(relative_change),
    max_abs_relative_change = max(abs(relative_change)),
    .groups = "drop"
  )
write.csv(smearing_summary,
          file.path(RESULTS_DIR, "simulation_smearing_summary.csv"), row.names = FALSE)

cat("\n=== Convergence summary ===\n")
print(as.data.frame(convergence_summary), row.names = FALSE, digits = 4)
cat("\n=== Smearing homoskedasticity summary ===\n")
print(as.data.frame(smearing_summary), row.names = FALSE, digits = 4)
cat("\nDiagnostics completed.\n")
