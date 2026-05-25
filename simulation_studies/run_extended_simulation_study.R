################################################################################
##  Extended Standard BCF Simulation Study for Semicontinuous Data
##
##  This script evaluates standard Gaussian BCF across three semicontinuous DGPs:
##    1. DGP A: Log-Normal Hurdle DGP (proposal standard)
##    2. DGP B: Gaussian (Normal) Hurdle DGP (continuous component is normal)
##    3. DGP C: Intrinsic Zero-Inflated DGP (Tweedie compound Poisson-Gamma, p=1.5)
##
##  For each DGP, we fit and compare:
##    - BCF-Linear: Standard BCF fitted directly on the raw semicontinuous Y
##    - BCF-Log: Standard BCF fitted on log(Y + 1) and re-transformed
################################################################################

library(countbcf, lib.loc = "local_lib")

set.seed(42)

## ---- Global specs -----------------------------------------------------------
N     <- 1000
P     <- 5
NBURN <- 500
NSIM  <- 1000
NTHIN <- 1

# Generate standard normal covariates
X <- matrix(rnorm(N * P), N, P)
colnames(X) <- paste0("X", 1:P)

# 1. Propensity Score & Treatment Assignment (Confounded, shared by all DGPs)
pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
Z    <- rbinom(N, 1, pi_x)

RESULTS_DIR <- "simulation_studies/results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

## ---- DGP A: Log-Normal Hurdle DGP -------------------------------------------
generate_dgp_lognormal <- function() {
  p_hurdle_0   <- pnorm(0.2 + 0.5 * X[, 1] - 0.3 * X[, 3])
  p_hurdle_1   <- pnorm(0.2 + 0.5 * X[, 1] - 0.3 * X[, 3] + 0.4 + 0.2 * X[, 1])
  p_hurdle_obs <- ifelse(Z == 1, p_hurdle_1, p_hurdle_0)
  
  I <- rbinom(N, 1, p_hurdle_obs)
  
  mu_c_0     <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4]
  mu_c_1     <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4] + 0.5 - 0.3 * X[, 2]
  sigma_true <- 0.5
  
  y_pos_0   <- exp(mu_c_0 + rnorm(N, 0, sigma_true))
  y_pos_1   <- exp(mu_c_1 + rnorm(N, 0, sigma_true))
  y_pos_obs <- ifelse(Z == 1, y_pos_1, y_pos_0)
  
  Y <- I * y_pos_obs
  
  true_mu0  <- p_hurdle_0 * exp(mu_c_0 + 0.5 * sigma_true^2)
  true_mu1  <- p_hurdle_1 * exp(mu_c_1 + 0.5 * sigma_true^2)
  true_cate <- true_mu1 - true_mu0
  
  list(y = Y, true_cate = true_cate, true_ate = mean(true_cate), name = "Log-Normal Hurdle")
}

## ---- DGP B: Gaussian Hurdle DGP --------------------------------------------
generate_dgp_gaussian <- function() {
  p_hurdle_0   <- pnorm(0.2 + 0.5 * X[, 1] - 0.3 * X[, 3])
  p_hurdle_1   <- pnorm(0.2 + 0.5 * X[, 1] - 0.3 * X[, 3] + 0.4 + 0.2 * X[, 1])
  p_hurdle_obs <- ifelse(Z == 1, p_hurdle_1, p_hurdle_0)
  
  I <- rbinom(N, 1, p_hurdle_obs)
  
  # Higher mean to keep active values positive
  mu_c_0     <- 6.0 + 1.2 * X[, 2] + 0.6 * X[, 4]
  mu_c_1     <- 6.0 + 1.2 * X[, 2] + 0.6 * X[, 4] + 1.5 - 0.5 * X[, 2]
  sigma_true <- 1.0
  
  y_pos_0   <- mu_c_0 + rnorm(N, 0, sigma_true)
  y_pos_1   <- mu_c_1 + rnorm(N, 0, sigma_true)
  y_pos_obs <- ifelse(Z == 1, y_pos_1, y_pos_0)
  
  Y <- I * y_pos_obs
  
  true_mu0  <- p_hurdle_0 * mu_c_0
  true_mu1  <- p_hurdle_1 * mu_c_1
  true_cate <- true_mu1 - true_mu0
  
  list(y = Y, true_cate = true_cate, true_ate = mean(true_cate), name = "Gaussian Hurdle")
}

## ---- DGP C: Tweedie Compound Poisson-Gamma DGP (p=1.5) ----------------------
generate_dgp_tweedie <- function() {
  # True log-mean parameter
  log_mu0 <- 1.2 + 0.8 * X[, 1] - 0.4 * X[, 3]
  log_mu1 <- 1.2 + 0.8 * X[, 1] - 0.4 * X[, 3] + 0.6 + 0.3 * X[, 1]
  
  mu0_true  <- exp(log_mu0)
  mu1_true  <- exp(log_mu1)
  true_cate <- mu1_true - mu0_true
  
  mu_true  <- ifelse(Z == 1, mu1_true, mu0_true)
  phi_true <- 1.5
  
  # Simulation via latent Poisson and continuous Gamma sizes
  lambda_true <- 2 * sqrt(mu_true) / phi_true
  N_latent    <- rpois(N, lambda_true)
  gamma_true  <- 0.5 * phi_true * sqrt(mu_true)
  
  Y <- rep(0, N)
  for (i in 1:N) {
    if (N_latent[i] > 0) {
      Y[i] <- rgamma(1, shape = N_latent[i], scale = gamma_true[i])
    }
  }
  
  list(y = Y, true_cate = true_cate, true_ate = mean(true_cate), name = "Tweedie Compound")
}

## ---- Fitting & Evaluation Loop ----------------------------------------------
run_eval <- function(dgp) {
  cat(sprintf("\n=== Evaluating DGP: %s (True ATE = %.4f) ===\n", dgp$name, dgp$true_ate))
  
  # ----------------------------------------------------
  # Model 1: BCF-Linear (fitting on raw Y)
  # ----------------------------------------------------
  cat("  -> Fitting BCF-Linear...\n")
  fit_linear <- bcf_continuous_linear(
      y          = dgp$y,
      z          = Z,
      x_control  = X,
      x_moderate = X,
      zhat       = pi_x,
      nburn      = NBURN,
      nsim       = NSIM,
      nthin      = NTHIN,
      update_interval = 500
  )
  
  cate_draws_linear <- get_forest_fit(fit_linear$moderate_fit, X)
  ate_draws_linear  <- rowMeans(cate_draws_linear)
  
  cate_est_linear   <- colMeans(cate_draws_linear)
  cate_ci_linear    <- apply(cate_draws_linear, 2, quantile, probs = c(0.025, 0.975))
  
  rmse_linear     <- sqrt(mean((cate_est_linear - dgp$true_cate)^2))
  bias_linear     <- mean(cate_est_linear - dgp$true_cate)
  coverage_linear <- mean(dgp$true_cate >= cate_ci_linear[1, ] & dgp$true_cate <= cate_ci_linear[2, ])
  cor_linear      <- cor(cate_est_linear, dgp$true_cate)
  
  # ----------------------------------------------------
  # Model 2: BCF-Log (fitting on log(Y+1))
  # ----------------------------------------------------
  cat("  -> Fitting BCF-Log...\n")
  Y_log <- log(dgp$y + 1)
  muy_log <- mean(Y_log)
  fit_log <- bcf_continuous_linear(
      y          = Y_log,
      z          = Z,
      x_control  = X,
      x_moderate = X,
      zhat       = pi_x,
      nburn      = NBURN,
      nsim       = NSIM,
      nthin      = NTHIN,
      update_interval = 500
  )
  
  mu_post_log   <- muy_log + get_forest_fit(fit_log$control_fit, X)
  tau_post_log  <- get_forest_fit(fit_log$moderate_fit, X)
  sigma_post_log <- fit_log$sigma
  
  # Re-transform: Y = exp(W) - 1
  cate_draws_log <- matrix(0, nrow = NSIM, ncol = N)
  for (s in 1:NSIM) {
    mu0_draw <- exp(mu_post_log[s, ] + 0.5 * sigma_post_log[s]^2) - 1
    mu1_draw <- exp(mu_post_log[s, ] + tau_post_log[s, ] + 0.5 * sigma_post_log[s]^2) - 1
    cate_draws_log[s, ] <- mu1_draw - mu0_draw
  }
  ate_draws_log <- rowMeans(cate_draws_log)
  
  cate_est_log  <- colMeans(cate_draws_log)
  cate_ci_log   <- apply(cate_draws_log, 2, quantile, probs = c(0.025, 0.975))
  
  rmse_log     <- sqrt(mean((cate_est_log - dgp$true_cate)^2))
  bias_log     <- mean(cate_est_log - dgp$true_cate)
  coverage_log <- mean(dgp$true_cate >= cate_ci_log[1, ] & dgp$true_cate <= cate_ci_log[2, ])
  cor_log      <- cor(cate_est_log, dgp$true_cate)
  
  # Pack results
  data.frame(
    DGP = dgp$name,
    True_ATE = dgp$true_ate,
    Metric = c("Est ATE Mean", "Est ATE SD", "CATE RMSE", "CATE Abs Bias", "CATE 95% Coverage", "CATE Correlation"),
    BCF_Linear = c(mean(ate_draws_linear), sd(ate_draws_linear), rmse_linear, abs(bias_linear), coverage_linear, cor_linear),
    BCF_Log = c(mean(ate_draws_log), sd(ate_draws_log), rmse_log, abs(bias_log), coverage_log, cor_log)
  )
}

# Run all 3 DGPs
dgp_a <- generate_dgp_lognormal()
dgp_b <- generate_dgp_gaussian()
dgp_c <- generate_dgp_tweedie()

res_a <- run_eval(dgp_a)
res_b <- run_eval(dgp_b)
res_c <- run_eval(dgp_c)

# Combine and save results
all_results <- rbind(res_a, res_b, res_c)
write.csv(all_results, file.path(RESULTS_DIR, "extended_bcf_simulation_results.csv"), row.names = FALSE)
cat(sprintf("\n[SUCCESS] Extended simulation results saved to: %s\n", file.path(RESULTS_DIR, "extended_bcf_simulation_results.csv")))

# Print results
print(all_results, digits = 4)
