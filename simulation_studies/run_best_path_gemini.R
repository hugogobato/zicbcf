################################################################################
##  ZIC-BCF-Smear (Best_Path_Gemini) Simulation Study
##
##  This script evaluates the new ZIC-BCF-Smear model across three DGPs:
##    - DGP A: Log-Normal Hurdle DGP (proposal standard)
##    - DGP B: Gaussian (Normal) Hurdle DGP (continuous component is normal)
##    - DGP C: Intrinsic Zero-Inflated DGP (Tweedie compound Poisson-Gamma, p=1.5)
##
##  Evaluations are done under the same specs as the walkthrough:
##    - DGP A: nburn=500, nsim=1000
##    - DGP B: nburn=500, nsim=1000
##    - DGP C: nburn=500, nsim=1000
##    - DGP C Corrected: nburn=2000, nsim=1000
################################################################################

library(countbcf, lib.loc = "local_lib")

set.seed(42)

## ---- Global specs -----------------------------------------------------------
N     <- 1000
P     <- 5
NSIM  <- 1000
NTHIN <- 1

# Generate standard normal covariates
X <- matrix(rnorm(N * P), N, P)
colnames(X) <- paste0("X", 1:P)

# Propensity Score & Treatment Assignment (Confounded, shared by all DGPs)
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

## ---- Fitting & Evaluation Function ------------------------------------------
evaluate_smear <- function(dgp, nburn) {
  cat(sprintf("\n  -> Fitting ZIC-BCF-Smear on %s (nburn=%d)...\n", dgp$name, nburn))
  
  fit_smear <- zicbcf_smear(
      y             = dgp$y,
      z             = Z,
      x_control     = X,
      pihat         = pi_x,
      pihat_active  = NULL, # SPA
      nburn         = nburn,
      nsim          = NSIM,
      nthin         = NTHIN,
      update_interval = 9999
  )
  
  cate_draws_smear <- fit_smear$cate
  ate_draws_smear  <- fit_smear$ate
  
  cate_est_smear   <- colMeans(cate_draws_smear)
  cate_ci_smear    <- apply(cate_draws_smear, 2, quantile, probs = c(0.025, 0.975))
  
  rmse_smear     <- sqrt(mean((cate_est_smear - dgp$true_cate)^2))
  bias_smear     <- mean(cate_est_smear - dgp$true_cate)
  coverage_smear <- mean(dgp$true_cate >= cate_ci_smear[1, ] & dgp$true_cate <= cate_ci_smear[2, ])
  cor_smear      <- cor(cate_est_smear, dgp$true_cate)
  
  res_metrics <- c(
    mean(ate_draws_smear),
    sd(ate_draws_smear),
    rmse_smear,
    abs(bias_smear),
    coverage_smear,
    cor_smear
  )
  
  cat(sprintf("     Est ATE Mean = %.4f | CATE RMSE = %.4f | CATE Abs Bias = %.4f | Coverage = %.1f%% | Correlation = %.4f\n",
              res_metrics[1], res_metrics[3], res_metrics[4], res_metrics[5]*100, res_metrics[6]))
              
  return(res_metrics)
}

## ---- Main Simulation --------------------------------------------------------
cat("=== Running ZIC-BCF-Smear Simulation Studies ===\n")

dgp_a <- generate_dgp_lognormal()
dgp_b <- generate_dgp_gaussian()
dgp_c <- generate_dgp_tweedie()

res_dgp_a      <- evaluate_smear(dgp_a, nburn = 500)
res_dgp_b      <- evaluate_smear(dgp_b, nburn = 500)
res_dgp_c      <- evaluate_smear(dgp_c, nburn = 500)
res_dgp_c_2000 <- evaluate_smear(dgp_c, nburn = 2000)

out_df <- data.frame(
  Metric = c("Est ATE Mean", "Est ATE SD", "CATE RMSE", "CATE Abs Bias", "CATE 95% Coverage", "CATE Correlation"),
  DGP_A_nburn500 = res_dgp_a,
  DGP_B_nburn500 = res_dgp_b,
  DGP_C_nburn500 = res_dgp_c,
  DGP_C_nburn2000 = res_dgp_c_2000
)

print(out_df, digits = 4)
write.csv(out_df, "simulation_studies/results/best_path_gemini_results.csv", row.names = FALSE)
cat("\n[SUCCESS] ZIC-BCF-Smear simulation results saved to simulation_studies/results/best_path_gemini_results.csv\n")
