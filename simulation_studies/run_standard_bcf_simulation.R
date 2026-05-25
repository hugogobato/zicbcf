################################################################################
##  Standard BCF Simulation Study for Semicontinuous (Zero-Inflated Continuous) Data
##
##  This script implements a simulation study to test how standard Gaussian BCF
##  performs on zero-inflated continuous (semicontinuous) data under two settings:
##    1. Raw Linear Scale: fitting BCF directly on Y (which has a spike at 0).
##    2. Log-transformed Scale: fitting BCF on log(Y + 1) and re-transforming.
##
##  The Data-Generating Process (DGP) matches Section 4 of the ZIC-BCF proposal.
################################################################################

library(countbcf, lib.loc = "local_lib")

set.seed(42)

## ---- Simulation settings ----------------------------------------------------
N     <- 1000
P     <- 5
NBURN <- 500
NSIM  <- 1000
NTHIN <- 1

cat("========================================================================\n")
cat("1. Generating Semicontinuous Data according to Proposal DGP...\n")
cat("========================================================================\n")

# Generate standard normal covariates
X <- matrix(rnorm(N * P), N, P)
colnames(X) <- paste0("X", 1:P)

# 1. Propensity Score & Treatment Assignment (Confounded)
pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
Z    <- rbinom(N, 1, pi_x)

# 2. Binary Participation Hurdle (true hurdle probabilities)
p_hurdle_0   <- pnorm(0.2 + 0.5 * X[, 1] - 0.3 * X[, 3])
p_hurdle_1   <- pnorm(0.2 + 0.5 * X[, 1] - 0.3 * X[, 3] + 0.4 + 0.2 * X[, 1])
p_hurdle_obs <- ifelse(Z == 1, p_hurdle_1, p_hurdle_0)

I <- rbinom(N, 1, p_hurdle_obs)

# 3. Continuous Intensity Component (Log-Normal)
mu_c_0     <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4]
mu_c_1     <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4] + 0.5 - 0.3 * X[, 2]
sigma_true <- 0.5 # sd of error is sqrt(0.25) = 0.5

y_pos_0   <- exp(mu_c_0 + rnorm(N, 0, sigma_true))
y_pos_1   <- exp(mu_c_1 + rnorm(N, 0, sigma_true))
y_pos_obs <- ifelse(Z == 1, y_pos_1, y_pos_0)

# 4. Final Semicontinuous Outcome (Zero-Inflated Continuous)
Y <- I * y_pos_obs

# Compute true potential outcomes and CATE on the original response scale
true_mu0  <- p_hurdle_0 * exp(mu_c_0 + 0.5 * sigma_true^2)
true_mu1  <- p_hurdle_1 * exp(mu_c_1 + 0.5 * sigma_true^2)
true_cate <- true_mu1 - true_mu0
true_ate  <- mean(true_cate)

cat("Data Summary:\n")
cat(sprintf("  - Total observations:             %d\n", N))
cat(sprintf("  - Active (positive) proportion:  %.1f%%\n", 100 * mean(Y > 0)))
cat(sprintf("  - Treatment proportion (z = 1):  %.1f%%\n", 100 * mean(Z)))
cat(sprintf("  - Outcome Mean (Overall):         %.2f\n", mean(Y)))
cat(sprintf("  - Outcome Max (Overall):          %.2f\n", max(Y)))
cat(sprintf("  - True ATE (Response Scale):      %.4f\n", true_ate))
cat("\n")

## ---- Fit Method 1: Standard BCF on Raw Y ------------------------------------
cat("========================================================================\n")
cat("2. Fitting Standard BCF on Raw Y (Linear Scale)...\n")
cat("========================================================================\n")

t0_linear <- Sys.time()
fit_linear <- bcf_continuous_linear(
    y          = Y,
    z          = Z,
    x_control  = X,
    x_moderate = X,
    zhat       = pi_x,
    nburn      = NBURN,
    nsim       = NSIM,
    nthin      = NTHIN,
    update_interval = 250
)
t1_linear <- Sys.time()
elapsed_linear <- as.numeric(difftime(t1_linear, t0_linear, units = "secs"))
cat(sprintf("  Linear BCF complete in %.1f seconds.\n\n", elapsed_linear))

# Recover linear CATE draws (get_forest_fit automatically multiplies by sdy)
cate_draws_linear <- get_forest_fit(fit_linear$moderate_fit, X) # nsim x n
ate_draws_linear  <- rowMeans(cate_draws_linear)

# Summaries
cate_est_linear <- colMeans(cate_draws_linear)
cate_ci_linear  <- apply(cate_draws_linear, 2, quantile, probs = c(0.025, 0.975))

rmse_linear     <- sqrt(mean((cate_est_linear - true_cate)^2))
bias_linear     <- mean(cate_est_linear - true_cate)
coverage_linear <- mean(true_cate >= cate_ci_linear[1, ] & true_cate <= cate_ci_linear[2, ])
cor_linear      <- cor(cate_est_linear, true_cate)

## ---- Fit Method 2: Standard BCF on log(Y + 1) --------------------------------
cat("========================================================================\n")
cat("3. Fitting Standard BCF on log(Y + 1) (Log Scale)...\n")
cat("========================================================================\n")

Y_log <- log(Y + 1)
muy_log <- mean(Y_log)

t0_log <- Sys.time()
fit_log <- bcf_continuous_linear(
    y          = Y_log,
    z          = Z,
    x_control  = X,
    x_moderate = X,
    zhat       = pi_x,
    nburn      = NBURN,
    nsim       = NSIM,
    nthin      = NTHIN,
    update_interval = 250
)
t1_log <- Sys.time()
elapsed_log <- as.numeric(difftime(t1_log, t0_log, units = "secs"))
cat(sprintf("  Log BCF complete in %.1f seconds.\n\n", elapsed_log))

# Recover log-scale posteriors
mu_post_log   <- muy_log + get_forest_fit(fit_log$control_fit, X)   # nsim x n
tau_post_log  <- get_forest_fit(fit_log$moderate_fit, X)            # nsim x n
sigma_post_log <- fit_log$sigma                                      # nsim

# Re-transform draws to the response scale: Y = exp(W) - 1
cate_draws_log <- matrix(0, nrow = NSIM, ncol = N)
for (s in 1:NSIM) {
  mu0_draw <- exp(mu_post_log[s, ] + 0.5 * sigma_post_log[s]^2) - 1
  mu1_draw <- exp(mu_post_log[s, ] + tau_post_log[s, ] + 0.5 * sigma_post_log[s]^2) - 1
  cate_draws_log[s, ] <- mu1_draw - mu0_draw
}
ate_draws_log <- rowMeans(cate_draws_log)

# Summaries
cate_est_log <- colMeans(cate_draws_log)
cate_ci_log  <- apply(cate_draws_log, 2, quantile, probs = c(0.025, 0.975))

rmse_log     <- sqrt(mean((cate_est_log - true_cate)^2))
bias_log     <- mean(cate_est_log - true_cate)
coverage_log <- mean(true_cate >= cate_ci_log[1, ] & true_cate <= cate_ci_log[2, ])
cor_log      <- cor(cate_est_log, true_cate)

## ---- Print and Save Results -------------------------------------------------
cat("========================================================================\n")
cat("4. COMPARATIVE EVALUATION RESULTS (RESPONSE SCALE CATE)\n")
cat("========================================================================\n")

cat(sprintf("True ATE:                    %.4f\n\n", true_ate))

results_summary <- data.frame(
    Metric = c("Estimated ATE (Mean)", "Estimated ATE (SD)", "ATE 95% Credible Interval", 
               "CATE RMSE", "CATE Absolute Bias", "CATE 95% Coverage Rate", "CATE correlation"),
    BCF_Linear = c(
        sprintf("%.4f", mean(ate_draws_linear)),
        sprintf("%.4f", sd(ate_draws_linear)),
        sprintf("[%.4f, %.4f]", quantile(ate_draws_linear, 0.025), quantile(ate_draws_linear, 0.975)),
        sprintf("%.4f", rmse_linear),
        sprintf("%.4f", abs(bias_linear)),
        sprintf("%.1f%%", 100 * coverage_linear),
        sprintf("%.4f", cor_linear)
    ),
    BCF_Log = c(
        sprintf("%.4f", mean(ate_draws_log)),
        sprintf("%.4f", sd(ate_draws_log)),
        sprintf("[%.4f, %.4f]", quantile(ate_draws_log, 0.025), quantile(ate_draws_log, 0.975)),
        sprintf("%.4f", rmse_log),
        sprintf("%.4f", abs(bias_log)),
        sprintf("%.1f%%", 100 * coverage_log),
        sprintf("%.4f", cor_log)
    ),
    stringsAsFactors = FALSE
)

print(results_summary, row.names = FALSE)
cat("========================================================================\n")

# Save results
RESULTS_DIR <- "example/results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
write.csv(results_summary, file.path(RESULTS_DIR, "standard_bcf_simulation_results.csv"), row.names = FALSE)
cat(sprintf("Saved results to: %s\n", file.path(RESULTS_DIR, "standard_bcf_simulation_results.csv")))

# Also save sample predictions
predictions_df <- data.frame(
    Unit = 1:N,
    Z = Z,
    Y = Y,
    True_CATE = true_cate,
    Linear_CATE_Est = cate_est_linear,
    Log_CATE_Est = cate_est_log
)
write.csv(predictions_df, file.path(RESULTS_DIR, "standard_bcf_unit_predictions.csv"), row.names = FALSE)
cat(sprintf("Saved unit predictions to: %s\n", file.path(RESULTS_DIR, "standard_bcf_unit_predictions.csv")))
cat("========================================================================\n")
