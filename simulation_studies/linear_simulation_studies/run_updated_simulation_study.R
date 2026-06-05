################################################################################
##  Updated BCF vs ZIC-BCF-Smear Semicontinuous Simulation Study (Local)
##
##  This script runs the local comparative simulation study evaluating:
##    1. BCF-Linear: Standard BCF fitted directly on the raw semicontinuous Y.
##    2. ZIC-BCF-Smear: Decoupled two-part hurdle BCF with Duan's Smearing.
##
##  Across three zero-inflated and right-skewed DGPs:
##    - DGP A: Log-Normal Hurdle DGP (halved treatment effect)
##    - DGP B: Gamma Hurdle DGP (new replacement right-skewed DGP, shape = 2.0)
##    - DGP C: Tweedie Compound Poisson-Gamma DGP (halved treatment effect)
################################################################################

# Robust package loading
if (!require("zicbcf", quietly = TRUE)) {
  library(countbcf, lib.loc = "local_lib")
}
library(dplyr)
library(tidyr)

# Create results directory if it does not exist
RESULTS_DIR <- "simulation_studies/results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

## ---- Global Simulation Parameters -------------------------------------------
N_SIM <- 100L  # 100 seeds for high statistical precision
N     <- 500L  # Standard sample size
P     <- 5L
NBURN <- 1000L # Convergence parameters
NSIM  <- 1000L
NTHIN <- 1L

# Generate standard normal covariates
set.seed(42)
X <- matrix(rnorm(N * P), N, P)
colnames(X) <- paste0("X", 1:P)

# Confounded Propensity Score & Treatment Assignment
pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
Z    <- rbinom(N, 1, pi_x)

## ---- DGP Generators with Adjusted (Halved) Treatment Effects ----------------

generate_dgp_a <- function(seed) {
  set.seed(seed * 1000 + 42)
  X_s <- matrix(rnorm(N * P), N, P)
  colnames(X_s) <- paste0("X", 1:P)
  pi_s <- pnorm(-0.5 + 0.4 * X_s[, 1] + 0.3 * X_s[, 2]^2)
  Z_s  <- rbinom(N, 1, pi_s)
  
  p_hurdle_0   <- pnorm(0.2 + 0.5 * X_s[, 1] - 0.3 * X_s[, 3])
  p_hurdle_1   <- pnorm(0.2 + 0.5 * X_s[, 1] - 0.3 * X_s[, 3] + 0.2 + 0.1 * X_s[, 1])
  p_hurdle_obs <- ifelse(Z_s == 1, p_hurdle_1, p_hurdle_0)
  I <- rbinom(N, 1, p_hurdle_obs)
  
  mu_c_0     <- 1.5 + 0.8 * X_s[, 2] + 0.4 * X_s[, 4]
  mu_c_1     <- 1.5 + 0.8 * X_s[, 2] + 0.4 * X_s[, 4] + 0.25 - 0.15 * X_s[, 2]
  sigma_true <- 0.5
  
  y_pos_0   <- exp(mu_c_0 + rnorm(N, 0, sigma_true))
  y_pos_1   <- exp(mu_c_1 + rnorm(N, 0, sigma_true))
  y_pos_obs <- ifelse(Z_s == 1, y_pos_1, y_pos_0)
  Y <- I * y_pos_obs
  
  true_mu0  <- p_hurdle_0 * exp(mu_c_0 + 0.5 * sigma_true^2)
  true_mu1  <- p_hurdle_1 * exp(mu_c_1 + 0.5 * sigma_true^2)
  true_cate <- true_mu1 - true_mu0
  true_hurdle_cate <- p_hurdle_1 - p_hurdle_0
  
  list(y = Y, z = Z_s, x = X_s, pihat = pi_s, true_cate = true_cate, true_ate = mean(true_cate), 
       true_hurdle_cate = true_hurdle_cate, true_hurdle_ate = mean(true_hurdle_cate))
}

generate_dgp_b <- function(seed) {
  set.seed(seed * 1000 + 42)
  X_s <- matrix(rnorm(N * P), N, P)
  colnames(X_s) <- paste0("X", 1:P)
  pi_s <- pnorm(-0.5 + 0.4 * X_s[, 1] + 0.3 * X_s[, 2]^2)
  Z_s  <- rbinom(N, 1, pi_s)
  
  p_hurdle_0   <- pnorm(0.2 + 0.5 * X_s[, 1] - 0.3 * X_s[, 3])
  p_hurdle_1   <- pnorm(0.2 + 0.5 * X_s[, 1] - 0.3 * X_s[, 3] + 0.2 + 0.1 * X_s[, 1])
  p_hurdle_obs <- ifelse(Z_s == 1, p_hurdle_1, p_hurdle_0)
  I <- rbinom(N, 1, p_hurdle_obs)
  
  log_mu_c_0 <- 1.5 + 0.8 * X_s[, 2] + 0.4 * X_s[, 4]
  log_mu_c_1 <- 1.5 + 0.8 * X_s[, 2] + 0.4 * X_s[, 4] + 0.25 - 0.15 * X_s[, 2]
  
  mu_c_0 <- exp(log_mu_c_0)
  mu_c_1 <- exp(log_mu_c_1)
  
  alpha <- 2.0
  scale_0 <- mu_c_0 / alpha
  scale_1 <- mu_c_1 / alpha
  
  y_pos_0 <- rgamma(N, shape = alpha, scale = scale_0)
  y_pos_1 <- rgamma(N, shape = alpha, scale = scale_1)
  y_pos_obs <- ifelse(Z_s == 1, y_pos_1, y_pos_0)
  Y <- I * y_pos_obs
  
  true_mu0  <- p_hurdle_0 * mu_c_0
  true_mu1  <- p_hurdle_1 * mu_c_1
  true_cate <- true_mu1 - true_mu0
  true_hurdle_cate <- p_hurdle_1 - p_hurdle_0
  
  list(y = Y, z = Z_s, x = X_s, pihat = pi_s, true_cate = true_cate, true_ate = mean(true_cate), 
       true_hurdle_cate = true_hurdle_cate, true_hurdle_ate = mean(true_hurdle_cate))
}

generate_dgp_c <- function(seed) {
  set.seed(seed * 1000 + 42)
  X_s <- matrix(rnorm(N * P), N, P)
  colnames(X_s) <- paste0("X", 1:P)
  pi_s <- pnorm(-0.5 + 0.4 * X_s[, 1] + 0.3 * X_s[, 2]^2)
  Z_s  <- rbinom(N, 1, pi_s)
  
  log_mu0 <- 1.2 + 0.8 * X_s[, 1] - 0.4 * X_s[, 3]
  log_mu1 <- 1.2 + 0.8 * X_s[, 1] - 0.4 * X_s[, 3] + 0.3 + 0.15 * X_s[, 1]
  
  mu0_true  <- exp(log_mu0)
  mu1_true  <- exp(log_mu1)
  true_cate <- mu1_true - mu0_true
  
  mu_true  <- ifelse(Z_s == 1, mu1_true, mu0_true)
  phi_true <- 1.5
  
  lambda0_true <- 2 * sqrt(mu0_true) / phi_true
  lambda1_true <- 2 * sqrt(mu1_true) / phi_true
  
  lambda_true <- ifelse(Z_s == 1, lambda1_true, lambda0_true)
  N_latent    <- rpois(N, lambda_true)
  gamma_true  <- 0.5 * phi_true * sqrt(mu_true)
  
  Y <- rep(0, N)
  for (i in 1:N) {
    if (N_latent[i] > 0) {
      Y[i] <- rgamma(1, shape = N_latent[i], scale = gamma_true[i])
    }
  }
  
  p0_hurdle_true <- 1 - exp(-lambda0_true)
  p1_hurdle_true <- 1 - exp(-lambda1_true)
  true_hurdle_cate <- p1_hurdle_true - p0_hurdle_true
  
  list(y = Y, z = Z_s, x = X_s, pihat = pi_s, true_cate = true_cate, true_ate = mean(true_cate), 
       true_hurdle_cate = true_hurdle_cate, true_hurdle_ate = mean(true_hurdle_cate))
}

# Metric calculation helper
calc_cate_metrics <- function(cate_draws, true_c, ate_draws) {
  cate_est <- colMeans(cate_draws)
  cate_ci  <- apply(cate_draws, 2, quantile, probs = c(0.025, 0.975))
  
  rmse <- sqrt(mean((cate_est - true_c)^2))
  bias <- mean(cate_est - true_c)
  coverage <- mean(true_c >= cate_ci[1, ] & true_c <= cate_ci[2, ])
  correlation <- cor(cate_est, true_c)
  if (is.na(correlation)) correlation <- 0.0
  ci_length <- mean(cate_ci[2, ] - cate_ci[1, ])
  est_ate_mean <- mean(ate_draws)
  
  list(rmse=rmse, bias=bias, coverage=coverage, correlation=correlation, ci_length=ci_length, est_ate_mean=est_ate_mean)
}

## ---- Main Simulation Loop ---------------------------------------------------

dgp_list <- list(
  list(name = "DGP A: Log-Normal Hurdle", gen = generate_dgp_a, code = "A"),
  list(name = "DGP B: Gamma Hurdle", gen = generate_dgp_b, code = "B"),
  list(name = "DGP C: Tweedie Semicontinuous", gen = generate_dgp_c, code = "C")
)

results_accum <- list()
idx <- 1

OUT_CSV <- file.path(RESULTS_DIR, "updated_simulation_results.csv")
if (file.exists(OUT_CSV)) file.remove(OUT_CSV)

for (dgp_info in dgp_list) {
  cat(sprintf("\n========================================================================\n"))
  cat(sprintf("=== Running Study for %s ===\n", dgp_info$name))
  cat(sprintf("========================================================================\n"))
  
  for (s in 1:N_SIM) {
    cat(sprintf("  [Seed %d/%d] Generating and fitting...\\n", s, N_SIM))
    d <- dgp_info$gen(s)
    
    # 1. Fit BCF-Linear
    fit_lin <- bcf_continuous_linear(
      y          = d$y,
      z          = d$z,
      x_control  = d$x,
      x_moderate = d$x,
      zhat       = d$pihat,
      nburn      = NBURN,
      nsim       = NSIM,
      thin       = NTHIN,
      update_interval = 99999
    )
    cate_draws_lin <- get_forest_fit(fit_lin$moderate_fit, d$x)
    ate_draws_lin  <- rowMeans(cate_draws_lin)
    m_lin <- calc_cate_metrics(cate_draws_lin, d$true_cate, ate_draws_lin)
    
    # 2. Fit ZIC-BCF-Smear
    fit_smear <- zicbcf_smear(
      y             = d$y,
      z             = d$z,
      x_control     = d$x,
      x_moderate    = d$x,
      pihat         = d$pihat,
      nburn         = NBURN,
      nsim          = NSIM,
      update_interval = 99999
    )
    m_smear <- calc_cate_metrics(fit_smear$cate, d$true_cate, fit_smear$ate)
    
    # Probit Hurdle stage draws
    p0_draws <- pnorm(fit_smear$mu_b)
    p1_draws <- pnorm(fit_smear$mu_b + fit_smear$tau_b)
    hurdle_cate_draws <- p1_draws - p0_draws
    hurdle_ate_draws  <- rowMeans(hurdle_cate_draws)
    m_hurdle_smear <- calc_cate_metrics(hurdle_cate_draws, d$true_hurdle_cate, hurdle_ate_draws)
    
    df_row <- data.frame(
      DGP = dgp_info$name,
      Seed = s,
      True_ATE = d$true_ate,
      True_Hurdle_ATE = d$true_hurdle_ate,
      
      # Standard CATE Metrics
      Linear_CATE_RMSE        = m_lin$rmse,
      Linear_CATE_Abs_Bias    = abs(m_lin$bias),
      Linear_CATE_Coverage    = m_lin$coverage,
      Linear_CATE_Correlation = m_lin$correlation,
      Linear_CATE_CI_Length   = m_lin$ci_length,
      Linear_Est_ATE          = m_lin$est_ate_mean,
      
      PathA_CATE_RMSE         = NA,
      PathA_CATE_Abs_Bias     = NA,
      PathA_CATE_Coverage     = NA,
      PathA_CATE_Correlation  = NA,
      PathA_CATE_CI_Length    = NA,
      PathA_Est_ATE           = NA,
      
      Smear_CATE_RMSE        = m_smear$rmse,
      Smear_CATE_Abs_Bias    = abs(m_smear$bias),
      Smear_CATE_Coverage    = m_smear$coverage,
      Smear_CATE_Correlation = m_smear$correlation,
      Smear_CATE_CI_Length   = m_smear$ci_length,
      Smear_Est_ATE          = m_smear$est_ate_mean,
      
      Linear_Hurdle_RMSE        = NA,
      Linear_Hurdle_Abs_Bias    = NA,
      Linear_Hurdle_Coverage    = NA,
      Linear_Hurdle_Correlation = NA,
      Linear_Hurdle_CI_Length   = NA,
      Linear_Est_Hurdle_ATE     = NA,
      
      PathA_Hurdle_RMSE        = NA,
      PathA_Hurdle_Abs_Bias    = NA,
      PathA_Hurdle_Coverage    = NA,
      PathA_Hurdle_Correlation = NA,
      PathA_Hurdle_CI_Length   = NA,
      PathA_Est_Hurdle_ATE     = NA,
      
      Smear_Hurdle_RMSE        = m_hurdle_smear$rmse,
      Smear_Hurdle_Abs_Bias    = abs(m_hurdle_smear$bias),
      Smear_Hurdle_Coverage    = m_hurdle_smear$coverage,
      Smear_Hurdle_Correlation = m_hurdle_smear$correlation,
      Smear_Hurdle_CI_Length   = m_hurdle_smear$ci_length,
      Smear_Est_Hurdle_ATE     = m_hurdle_smear$est_ate_mean,
      
      stringsAsFactors = FALSE
    )
    
    write.table(df_row, OUT_CSV, sep=",", row.names=FALSE, col.names=!file.exists(OUT_CSV), append=TRUE)
    gc(verbose = FALSE)
  }
}
cat(sprintf("\n[SUCCESS] Master results CSV saved to: %s\n", OUT_CSV))
