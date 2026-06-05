################################################################################
##  Updated BCF vs ZIC-BCF-Smear Sensitivity Analyses (Local)
##
##  This script runs the local sensitivity analyses matching the Colab setup:
##    1. Sensitivity to Sample Size N (100, 250, 500, 1000)
##    2. Sensitivity to Zero-Inflation Proportion (varying hurdle intercepts)
##
##  All runs use N_SIM = 100, NBURN = 1000, NSIM = 1000.
################################################################################

# Robust package loading
if (!require("zicbcf", quietly = TRUE)) {
  library(countbcf, lib.loc = "local_lib")
}
library(dplyr)
library(tidyr)

RESULTS_DIR <- "simulation_studies/results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

## ---- Global Simulation Parameters -------------------------------------------
N_SIM  <- 100L
P      <- 5L
NBURN  <- 1000L
NSIM   <- 1000L
NTHIN  <- 1L

# Generate standard normal covariates
set.seed(42)
X_base <- matrix(rnorm(2000 * P), 2000, P) # base matrix up to N=2000
colnames(X_base) <- paste0("X", 1:P)
pi_x_base <- pnorm(-0.5 + 0.4 * X_base[, 1] + 0.3 * X_base[, 2]^2)
Z_base    <- rbinom(2000, 1, pi_x_base)

## ---- DGP Generators taking N and c_shift ------------------------------------

dgp_a_gen <- function(n, seed, c_shift = 0.2) {
  set.seed(seed * 1000 + 42)
  X <- matrix(rnorm(n * P), n, P)
  colnames(X) <- paste0("X", 1:P)
  pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
  Z    <- rbinom(n, 1, pi_x)
  
  p_hurdle_0   <- pnorm(c_shift + 0.5 * X[, 1] - 0.3 * X[, 3])
  p_hurdle_1   <- pnorm(c_shift + 0.5 * X[, 1] - 0.3 * X[, 3] + 0.2 + 0.1 * X[, 1])
  p_hurdle_obs <- ifelse(Z == 1, p_hurdle_1, p_hurdle_0)
  I <- rbinom(n, 1, p_hurdle_obs)
  
  mu_c_0     <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4]
  mu_c_1     <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4] + 0.25 - 0.15 * X[, 2]
  sigma_true <- 0.5
  
  y_pos_0   <- exp(mu_c_0 + rnorm(n, 0, sigma_true))
  y_pos_1   <- exp(mu_c_1 + rnorm(n, 0, sigma_true))
  y_pos_obs <- ifelse(Z == 1, y_pos_1, y_pos_0)
  Y <- I * y_pos_obs
  
  true_mu0  <- p_hurdle_0 * exp(mu_c_0 + 0.5 * sigma_true^2)
  true_mu1  <- p_hurdle_1 * exp(mu_c_1 + 0.5 * sigma_true^2)
  true_cate <- true_mu1 - true_mu0
  true_hurdle_cate <- p_hurdle_1 - p_hurdle_0
  
  list(y = Y, z = Z, x = X, pihat = pi_x, true_cate = true_cate, true_ate = mean(true_cate), 
       true_hurdle_cate = true_hurdle_cate, true_hurdle_ate = mean(true_hurdle_cate))
}

dgp_b_gen <- function(n, seed, c_shift = 0.2) {
  set.seed(seed * 1000 + 42)
  X <- matrix(rnorm(n * P), n, P)
  colnames(X) <- paste0("X", 1:P)
  pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
  Z    <- rbinom(n, 1, pi_x)
  
  p_hurdle_0   <- pnorm(c_shift + 0.5 * X[, 1] - 0.3 * X[, 3])
  p_hurdle_1   <- pnorm(c_shift + 0.5 * X[, 1] - 0.3 * X[, 3] + 0.2 + 0.1 * X[, 1])
  p_hurdle_obs <- ifelse(Z == 1, p_hurdle_1, p_hurdle_0)
  I <- rbinom(n, 1, p_hurdle_obs)
  
  log_mu_c_0 <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4]
  log_mu_c_1 <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4] + 0.25 - 0.15 * X[, 2]
  
  mu_c_0 <- exp(log_mu_c_0)
  mu_c_1 <- exp(log_mu_c_1)
  
  alpha <- 2.0
  scale_0 <- mu_c_0 / alpha
  scale_1 <- mu_c_1 / alpha
  
  y_pos_0 <- rgamma(n, shape = alpha, scale = scale_0)
  y_pos_1 <- rgamma(n, shape = alpha, scale = scale_1)
  y_pos_obs <- ifelse(Z == 1, y_pos_1, y_pos_0)
  Y <- I * y_pos_obs
  
  true_mu0  <- p_hurdle_0 * mu_c_0
  true_mu1  <- p_hurdle_1 * mu_c_1
  true_cate <- true_mu1 - true_mu0
  true_hurdle_cate <- p_hurdle_1 - p_hurdle_0
  
  list(y = Y, z = Z, x = X, pihat = pi_x, true_cate = true_cate, true_ate = mean(true_cate), 
       true_hurdle_cate = true_hurdle_cate, true_hurdle_ate = mean(true_hurdle_cate))
}

dgp_c_gen <- function(n, seed, c_shift = 0.0) {
  set.seed(seed * 1000 + 42)
  X <- matrix(rnorm(n * P), n, P)
  colnames(X) <- paste0("X", 1:P)
  pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
  Z    <- rbinom(n, 1, pi_x)
  
  log_mu0 <- 1.2 + c_shift + 0.8 * X[, 1] - 0.4 * X[, 3]
  log_mu1 <- 1.2 + c_shift + 0.8 * X[, 1] - 0.4 * X[, 3] + 0.3 + 0.15 * X[, 1]
  
  mu0_true  <- exp(log_mu0)
  mu1_true  <- exp(log_mu1)
  true_cate <- mu1_true - mu0_true
  
  mu_true  <- ifelse(Z == 1, mu1_true, mu0_true)
  phi_true <- 1.5
  
  lambda0_true <- 2 * sqrt(mu0_true) / phi_true
  lambda1_true <- 2 * sqrt(mu1_true) / phi_true
  
  lambda_true <- ifelse(Z == 1, lambda1_true, lambda0_true)
  N_latent    <- rpois(n, lambda_true)
  gamma_true  <- 0.5 * phi_true * sqrt(mu_true)
  
  Y <- rep(0, n)
  for (i in 1:n) {
    if (N_latent[i] > 0) {
      Y[i] <- rgamma(1, shape = N_latent[i], scale = gamma_true[i])
    }
  }
  
  p0_hurdle_true <- 1 - exp(-lambda0_true)
  p1_hurdle_true <- 1 - exp(-lambda1_true)
  true_hurdle_cate <- p1_hurdle_true - p0_hurdle_true
  
  list(y = Y, z = Z, x = X, pihat = pi_x, true_cate = true_cate, true_ate = mean(true_cate), 
       true_hurdle_cate = true_hurdle_cate, true_hurdle_ate = mean(true_hurdle_cate))
}

# Helpers
fit_bcf <- function(d) {
  fit <- bcf_continuous_linear(
    y          = d$y,
    z          = d$z,
    x_control  = d$x,
    x_moderate = d$x,
    zhat       = d$pihat,
    nburn      = NBURN,
    nsim       = NSIM,
    nthin      = NTHIN,
    update_interval = 99999
  )
  cate_draws <- get_forest_fit(fit$moderate_fit, d$x)
  ate_draws  <- rowMeans(cate_draws)
  list(cate_draws = cate_draws, ate_draws = ate_draws)
}

fit_smear <- function(d) {
  fit <- zicbcf_smear(
    y             = d$y,
    z             = d$z,
    x_control     = d$x,
    x_moderate    = d$x,
    pihat         = d$pihat,
    nburn         = NBURN,
    nsim          = NSIM,
    nthin         = NTHIN,
    update_interval = 99999
  )
  list(cate_draws = fit$cate, ate_draws = fit$ate)
}

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

## =============================================================================
##  SENSITIVITY ANALYSIS 1: SAMPLE SIZE N
## =============================================================================

cat("=== Starting Local N-Sensitivity Analysis ===\n")
N_LIST <- c(100L, 250L, 500L, 1000L)

dgp_scenarios <- list(
  list(name = "DGP A: Log-Normal Hurdle", gen = dgp_a_gen, default_shift = 0.2),
  list(name = "DGP B: Gamma Hurdle", gen = dgp_b_gen, default_shift = 0.2),
  list(name = "DGP C: Tweedie Semicontinuous", gen = dgp_c_gen, default_shift = 0.0)
)

for (n_val in N_LIST) {
  OUT_CSV <- file.path(RESULTS_DIR, sprintf("local_sensitivity_n%d_results.csv", n_val))
  if (file.exists(OUT_CSV)) file.remove(OUT_CSV)
  
  for (scenario in dgp_scenarios) {
    dgp_name <- scenario$name
    dgp_gen  <- scenario$gen
    c_shift  <- scenario$default_shift
    
    for (model_name in c("BCF-Linear", "ZIC-BCF-Smear")) {
      fit_fn <- if (model_name == "BCF-Linear") fit_bcf else fit_smear
      
      for (s in 1:N_SIM) {
        cat(sprintf("[%s | N=%d | %s] Seed %d/%d...\\n", dgp_name, n_val, model_name, s, N_SIM))
        d <- dgp_gen(n_val, seed = s, c_shift = c_shift)
        
        fit <- tryCatch(fit_fn(d), error = function(e) NULL)
        
        if (!is.null(fit)) {
          m <- calc_cate_metrics(fit$cate_draws, d$true_cate, fit$ate_draws)
          
          # Compute hurdle metrics if smear model
          if (model_name == "ZIC-BCF-Smear") {
            # zicbcf_smear returns fit in fit_smear but we need hurdle draws
            # let's just fit and get them manually to perfectly match
            fit_raw <- zicbcf_smear(
              y             = d$y,
              z             = d$z,
              x_control     = d$x,
              pihat         = d$pihat,
              nburn         = NBURN,
              nsim          = NSIM,
              update_interval = 99999
            )
            p0_draws <- pnorm(fit_raw$mu_b)
            p1_draws <- pnorm(fit_raw$mu_b + fit_raw$tau_b)
            hurdle_cate_draws <- p1_draws - p0_draws
            hurdle_ate_draws  <- rowMeans(hurdle_cate_draws)
            m_hurdle <- calc_cate_metrics(hurdle_cate_draws, d$true_hurdle_cate, hurdle_ate_draws)
            
            df_row <- data.frame(
              DGP = dgp_name,
              Seed = s,
              True_ATE = d$true_ate,
              True_Hurdle_ATE = d$true_hurdle_ate,
              
              Linear_CATE_RMSE         = NA,
              Linear_CATE_Abs_Bias     = NA,
              Linear_CATE_Coverage     = NA,
              Linear_CATE_Correlation  = NA,
              Linear_CATE_CI_Length    = NA,
              Linear_Est_ATE           = NA,
              
              Smear_CATE_RMSE        = m$rmse,
              Smear_CATE_Abs_Bias    = abs(m$bias),
              Smear_CATE_Coverage    = m$coverage,
              Smear_CATE_Correlation = m$correlation,
              Smear_CATE_CI_Length   = m$ci_length,
              Smear_Est_ATE          = m$est_ate_mean,
              
              Linear_Hurdle_RMSE        = NA,
              Linear_Hurdle_Abs_Bias    = NA,
              Linear_Hurdle_Coverage    = NA,
              Linear_Hurdle_Correlation = NA,
              Linear_Hurdle_CI_Length   = NA,
              Linear_Est_Hurdle_ATE     = NA,
              
              Smear_Hurdle_RMSE        = m_hurdle$rmse,
              Smear_Hurdle_Abs_Bias    = abs(m_hurdle$bias),
              Smear_Hurdle_Coverage    = m_hurdle$coverage,
              Smear_Hurdle_Correlation = m_hurdle$correlation,
              Smear_Hurdle_CI_Length   = m_hurdle$ci_length,
              Smear_Est_Hurdle_ATE     = m_hurdle$est_ate_mean,
              stringsAsFactors = FALSE
            )
          } else {
            df_row <- data.frame(
              DGP = dgp_name,
              Seed = s,
              True_ATE = d$true_ate,
              True_Hurdle_ATE = d$true_hurdle_ate,
              
              Linear_CATE_RMSE        = m$rmse,
              Linear_CATE_Abs_Bias    = abs(m$bias),
              Linear_CATE_Coverage    = m$coverage,
              Linear_CATE_Correlation = m$correlation,
              Linear_CATE_CI_Length   = m$ci_length,
              Linear_Est_ATE          = m$est_ate_mean,
              
              Smear_CATE_RMSE         = NA,
              Smear_CATE_Abs_Bias     = NA,
              Smear_CATE_Coverage     = NA,
              Smear_CATE_Correlation  = NA,
              Smear_CATE_CI_Length    = NA,
              Smear_Est_ATE           = NA,
              
              Linear_Hurdle_RMSE        = NA,
              Linear_Hurdle_Abs_Bias    = NA,
              Linear_Hurdle_Coverage    = NA,
              Linear_Hurdle_Correlation = NA,
              Linear_Hurdle_CI_Length   = NA,
              Linear_Est_Hurdle_ATE     = NA,
              
              Smear_Hurdle_RMSE        = NA,
              Smear_Hurdle_Abs_Bias    = NA,
              Smear_Hurdle_Coverage    = NA,
              Smear_Hurdle_Correlation = NA,
              Smear_Hurdle_CI_Length   = NA,
              Smear_Est_Hurdle_ATE     = NA,
              stringsAsFactors = FALSE
            )
          }
          
          write.table(df_row, OUT_CSV, sep=\",\", row.names=FALSE, col.names=!file.exists(OUT_CSV), append=TRUE)
        }
        gc(verbose = FALSE)
      }
    }
  }
}
cat("[SUCCESS] N-Sensitivity local results saved.\n")

## =============================================================================
##  SENSITIVITY ANALYSIS 2: ZERO-INFLATION PROPORTION
## =============================================================================

cat("\n=== Starting Local ZI-Sensitivity Analysis ===\n")
N_FIXED <- 500L

zi_scenarios <- list(
  1 = list(a_shift = -1.5, b_shift = -1.5, c_shift = -3.5),
  2 = list(a_shift = -0.5, b_shift = -0.5, c_shift = -2.0),
  4 = list(a_shift = 1.0,  b_shift = 1.0,  c_shift = 0.0),
  5 = list(a_shift = 1.8,  b_shift = 1.8,  c_shift = 1.0)
)

for (lvl in names(zi_scenarios)) {
  shifts <- zi_scenarios[[lvl]]
  OUT_CSV <- file.path(RESULTS_DIR, sprintf("local_sensitivity_zi_lvl%s_results.csv", lvl))
  if (file.exists(OUT_CSV)) file.remove(OUT_CSV)
  
  dgp_zi_configs <- list(
    list(name = "DGP A: Log-Normal Hurdle", gen = dgp_a_gen, shift = shifts$a_shift),
    list(name = "DGP B: Gamma Hurdle", gen = dgp_b_gen, shift = shifts$b_shift),
    list(name = "DGP C: Tweedie Semicontinuous", gen = dgp_c_gen, shift = shifts$c_shift)
  )
  
  for (cfg in dgp_zi_configs) {
    for (model_name in c("BCF-Linear", "ZIC-BCF-Smear")) {
      fit_fn <- if (model_name == "BCF-Linear") fit_bcf else fit_smear
      
      for (s in 1:N_SIM) {
        cat(sprintf("[%s | Level %s | %s] Seed %d/%d...\\n", cfg$name, lvl, model_name, s, N_SIM))
        d <- cfg$gen(N_FIXED, seed = s, c_shift = cfg$shift)
        
        fit <- tryCatch(fit_fn(d), error = function(e) NULL)
        
        if (!is.null(fit)) {
          m <- calc_cate_metrics(fit$cate_draws, d$true_cate, fit$ate_draws)
          
          if (model_name == "ZIC-BCF-Smear") {
            fit_raw <- zicbcf_smear(
              y             = d$y,
              z             = d$z,
              x_control     = d$x,
              pihat         = d$pihat,
              nburn         = NBURN,
              nsim          = NSIM,
              update_interval = 99999
            )
            p0_draws <- pnorm(fit_raw$mu_b)
            p1_draws <- pnorm(fit_raw$mu_b + fit_raw$tau_b)
            hurdle_cate_draws <- p1_draws - p0_draws
            hurdle_ate_draws  <- rowMeans(hurdle_cate_draws)
            m_hurdle <- calc_cate_metrics(hurdle_cate_draws, d$true_hurdle_cate, hurdle_ate_draws)
            
            df_row <- data.frame(
              DGP = cfg$name,
              Seed = s,
              True_ATE = d$true_ate,
              True_Hurdle_ATE = d$true_hurdle_ate,
              
              Linear_CATE_RMSE         = NA,
              Linear_CATE_Abs_Bias     = NA,
              Linear_CATE_Coverage     = NA,
              Linear_CATE_Correlation  = NA,
              Linear_CATE_CI_Length    = NA,
              Linear_Est_ATE           = NA,
              
              Smear_CATE_RMSE        = m$rmse,
              Smear_CATE_Abs_Bias    = abs(m$bias),
              Smear_CATE_Coverage    = m$coverage,
              Smear_CATE_Correlation = m$correlation,
              Smear_CATE_CI_Length   = m$ci_length,
              Smear_Est_ATE          = m$est_ate_mean,
              
              Linear_Hurdle_RMSE        = NA,
              Linear_Hurdle_Abs_Bias    = NA,
              Linear_Hurdle_Coverage    = NA,
              Linear_Hurdle_Correlation = NA,
              Linear_Hurdle_CI_Length   = NA,
              Linear_Est_Hurdle_ATE     = NA,
              
              Smear_Hurdle_RMSE        = m_hurdle$rmse,
              Smear_Hurdle_Abs_Bias    = abs(m_hurdle$bias),
              Smear_Hurdle_Coverage    = m_hurdle$coverage,
              Smear_Hurdle_Correlation = m_hurdle$correlation,
              Smear_Hurdle_CI_Length   = m_hurdle$ci_length,
              Smear_Est_Hurdle_ATE     = m_hurdle$est_ate_mean,
              stringsAsFactors = FALSE
            )
          } else {
            df_row <- data.frame(
              DGP = cfg$name,
              Seed = s,
              True_ATE = d$true_ate,
              True_Hurdle_ATE = d$true_hurdle_ate,
              
              Linear_CATE_RMSE        = m$rmse,
              Linear_CATE_Abs_Bias    = abs(m$bias),
              Linear_CATE_Coverage    = m$coverage,
              Linear_CATE_Correlation = m$correlation,
              Linear_CATE_CI_Length   = m$ci_length,
              Linear_Est_ATE          = m$est_ate_mean,
              
              Smear_CATE_RMSE         = NA,
              Smear_CATE_Abs_Bias     = NA,
              Smear_CATE_Coverage     = NA,
              Smear_CATE_Correlation  = NA,
              Smear_CATE_CI_Length    = NA,
              Smear_Est_ATE           = NA,
              
              Linear_Hurdle_RMSE        = NA,
              Linear_Hurdle_Abs_Bias    = NA,
              Linear_Hurdle_Coverage    = NA,
              Linear_Hurdle_Correlation = NA,
              Linear_Hurdle_CI_Length   = NA,
              Linear_Est_Hurdle_ATE     = NA,
              
              Smear_Hurdle_RMSE        = NA,
              Smear_Hurdle_Abs_Bias    = NA,
              Smear_Hurdle_Coverage    = NA,
              Smear_Hurdle_Correlation = NA,
              Smear_Hurdle_CI_Length   = NA,
              Smear_Est_Hurdle_ATE     = NA,
              stringsAsFactors = FALSE
            )
          }
          
          write.table(df_row, OUT_CSV, sep=\",\", row.names=FALSE, col.names=!file.exists(OUT_CSV), append=TRUE)
        }
        gc(verbose = FALSE)
      }
    }
  }
}
cat("[SUCCESS] ZI-Sensitivity local results saved.\n")
cat("\n=== All Sensitivity Analyses Completed Locally ===\n")
