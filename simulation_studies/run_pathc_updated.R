################################################################################
##  Retesting Path E (Improved Selection Model with Active Outcome Fit)
##
##  This script evaluates:
##    - "Path E": Collapsed Joint Copula-BCF (active-subset outcome forest)
##
##  Using seed (42) and the exact same data-generating processes (DGPs).
##  Appends results to full_simulation_results.csv and dgpc_nburn2000_results.csv.
################################################################################

library(countbcf, lib.loc = "local_lib")

set.seed(42)

## ---- Global specs -----------------------------------------------------------
N     <- 1000
P     <- 5
NSIM  <- 1000
NTHIN <- 1

# Generate covariates and treatment assignment (identical to original)
X <- matrix(rnorm(N * P), N, P)
colnames(X) <- paste0("X", 1:P)
pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
Z    <- rbinom(N, 1, pi_x)

RESULTS_DIR <- "simulation_studies/results"

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
  log_mu0 <- 1.2 + 0.8 * X[, 1] - 0.4 * X[, 3]
  log_mu1 <- 1.2 + 0.8 * X[, 1] - 0.4 * X[, 3] + 0.6 + 0.3 * X[, 1]
  
  mu0_true  <- exp(log_mu0)
  mu1_true  <- exp(log_mu1)
  true_cate <- mu1_true - mu0_true
  
  mu_true  <- ifelse(Z == 1, mu1_true, mu0_true)
  phi_true <- 1.5
  
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

# Generate datasets sequentially under seed 42 (exactly matching original seed sequence)
dgp_a <- generate_dgp_lognormal()
dgp_b <- generate_dgp_gaussian()
dgp_c <- generate_dgp_tweedie()

## ---- Evaluation Function: Path E --------------------------------------------
run_pathe <- function(dgp, nburn) {
  cat(sprintf("\n--- Fitting Path E | DGP: %s (nburn = %d) ---\n", dgp$name, nburn))
  
  fit_pathE <- pathe_bcf(
      y          = dgp$y,
      z          = Z,
      x_control  = X,
      pihat_sel  = pi_x,
      nburn      = nburn,
      nsim       = NSIM,
      nthin      = NTHIN,
      update_interval = 500
  )
  
  cate_draws <- matrix(0, nrow = NSIM, ncol = N)
  for (s in 1:NSIM) {
    eta_b0 <- fit_pathE$sel_con_post[s, ]
    eta_b1 <- fit_pathE$sel_con_post[s, ] + fit_pathE$sel_tau_post[s, ]
    
    eta_c0 <- fit_pathE$out_con_post[s, ]
    eta_c1 <- fit_pathE$out_con_post[s, ] + fit_pathE$out_tau_post[s, ]
    
    sig2 <- fit_pathE$sigma_post[s]^2
    bet <- fit_pathE$beta_post[s]
    
    mu0_draw <- exp(eta_c0 + 0.5 * sig2) * pnorm(eta_b0 + bet)
    mu1_draw <- exp(eta_c1 + 0.5 * sig2) * pnorm(eta_b1 + bet)
    
    cate_draws[s, ] <- mu1_draw - mu0_draw
  }
  ate_draws <- rowMeans(cate_draws)
  cate_est  <- colMeans(cate_draws)
  cate_ci   <- apply(cate_draws, 2, quantile, probs = c(0.025, 0.975))
  
  rmse     <- sqrt(mean((cate_est - dgp$true_cate)^2))
  bias     <- mean(cate_est - dgp$true_cate)
  coverage <- mean(dgp$true_cate >= cate_ci[1, ] & dgp$true_cate <= cate_ci[2, ])
  cor_val  <- cor(cate_est, dgp$true_cate)
  
  list(
    ate_mean = mean(ate_draws),
    ate_sd = sd(ate_draws),
    rmse = rmse,
    bias = abs(bias),
    coverage = coverage,
    correlation = cor_val
  )
}

## ---- Run the evaluations ----------------------------------------------------
cat("\nRunning evaluations for Path E...\n")
res_e_a <- run_pathe(dgp_a, nburn = 500)
res_e_b <- run_pathe(dgp_b, nburn = 500)
res_e_c <- run_pathe(dgp_c, nburn = 500)
res_e_c_2000 <- run_pathe(dgp_c, nburn = 2000)

## ---- Append to full_simulation_results.csv ----------------------------------
csv_file_full <- file.path(RESULTS_DIR, "full_simulation_results.csv")
if (file.exists(csv_file_full)) {
  results_df <- read.csv(csv_file_full)
  
  # Path E column
  Path_E <- rep(0, nrow(results_df))
  # DGP A
  Path_E[results_df$DGP == "Log-Normal Hurdle" & results_df$Metric == "Est ATE Mean"] <- res_e_a$ate_mean
  Path_E[results_df$DGP == "Log-Normal Hurdle" & results_df$Metric == "Est ATE SD"]   <- res_e_a$ate_sd
  Path_E[results_df$DGP == "Log-Normal Hurdle" & results_df$Metric == "CATE RMSE"]    <- res_e_a$rmse
  Path_E[results_df$DGP == "Log-Normal Hurdle" & results_df$Metric == "CATE Abs Bias"]<- res_e_a$bias
  Path_E[results_df$DGP == "Log-Normal Hurdle" & results_df$Metric == "CATE 95% Coverage"] <- res_e_a$coverage
  Path_E[results_df$DGP == "Log-Normal Hurdle" & results_df$Metric == "CATE Correlation"]  <- res_e_a$correlation
  # DGP B
  Path_E[results_df$DGP == "Gaussian Hurdle" & results_df$Metric == "Est ATE Mean"] <- res_e_b$ate_mean
  Path_E[results_df$DGP == "Gaussian Hurdle" & results_df$Metric == "Est ATE SD"]   <- res_e_b$ate_sd
  Path_E[results_df$DGP == "Gaussian Hurdle" & results_df$Metric == "CATE RMSE"]    <- res_e_b$rmse
  Path_E[results_df$DGP == "Gaussian Hurdle" & results_df$Metric == "CATE Abs Bias"]<- res_e_b$bias
  Path_E[results_df$DGP == "Gaussian Hurdle" & results_df$Metric == "CATE 95% Coverage"] <- res_e_b$coverage
  Path_E[results_df$DGP == "Gaussian Hurdle" & results_df$Metric == "CATE Correlation"]  <- res_e_b$correlation
  # DGP C
  Path_E[results_df$DGP == "Tweedie Compound" & results_df$Metric == "Est ATE Mean"] <- res_e_c$ate_mean
  Path_E[results_df$DGP == "Tweedie Compound" & results_df$Metric == "Est ATE SD"]   <- res_e_c$ate_sd
  Path_E[results_df$DGP == "Tweedie Compound" & results_df$Metric == "CATE RMSE"]    <- res_e_c$rmse
  Path_E[results_df$DGP == "Tweedie Compound" & results_df$Metric == "CATE Abs Bias"]<- res_e_c$bias
  Path_E[results_df$DGP == "Tweedie Compound" & results_df$Metric == "CATE 95% Coverage"] <- res_e_c$coverage
  Path_E[results_df$DGP == "Tweedie Compound" & results_df$Metric == "CATE Correlation"]  <- res_e_c$correlation
  
  results_df$Path_E <- Path_E
  
  write.csv(results_df, csv_file_full, row.names = FALSE)
  cat("\n[SUCCESS] Updated full_simulation_results.csv with Path_E column.\n")
  print(results_df, digits = 4)
}

## ---- Append to dgpc_nburn2000_results.csv -----------------------------------
csv_file_2000 <- file.path(RESULTS_DIR, "dgpc_nburn2000_results.csv")
if (file.exists(csv_file_2000)) {
  results_2000_df <- read.csv(csv_file_2000)
  
  # Path E
  Path_E <- rep(0, nrow(results_2000_df))
  Path_E[results_2000_df$Metric == "Est ATE Mean"] <- res_e_c_2000$ate_mean
  Path_E[results_2000_df$Metric == "Est ATE SD"]   <- res_e_c_2000$ate_sd
  Path_E[results_2000_df$Metric == "CATE RMSE"]    <- res_e_c_2000$rmse
  Path_E[results_2000_df$Metric == "CATE Abs Bias"]<- res_e_c_2000$bias
  Path_E[results_2000_df$Metric == "CATE 95% Coverage"] <- res_e_c_2000$coverage
  Path_E[results_2000_df$Metric == "CATE Correlation"]  <- res_e_c_2000$correlation
  results_2000_df$Path_E <- Path_E
  
  write.csv(results_2000_df, csv_file_2000, row.names = FALSE)
  cat("\n[SUCCESS] Updated dgpc_nburn2000_results.csv with Path_E column.\n")
  print(results_2000_df, digits = 4)
}
