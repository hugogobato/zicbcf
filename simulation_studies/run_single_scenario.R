################################################################################
##  Zero-Inflation Sensitivity Analysis: Single Scenario Executor
##
##  This script runs a single scenario (1-15) in isolation to ensure that
##  the operating system reclaims 100% of the allocated memory upon completion.
##
##  Usage:
##    Rscript simulation_studies/run_single_scenario.R <scenario_id>
################################################################################

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Error: Must specify a scenario ID (1 to 15) as a command line argument.")
}
scenario_id <- as.integer(args[1])
if (is.na(scenario_id) || scenario_id < 1 || scenario_id > 15) {
  stop("Error: Scenario ID must be an integer between 1 and 15.")
}

library(countbcf, lib.loc = "local_lib")
library(tidyr)
library(dplyr)

# Create results directory
RESULTS_DIR <- "simulation_studies/results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

## ---- Global Simulation Specs -------------------------------------------------
N     <- 1000
P     <- 5
NBURN <- 500
NSIM  <- 1000
NTHIN <- 1

# Generate standard normal covariates
set.seed(42)
X <- matrix(rnorm(N * P), N, P)
colnames(X) <- paste0("X", 1:P)

# Propensity Score & Treatment Assignment (Confounded, shared by all runs)
pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
Z    <- rbinom(N, 1, pi_x)

## ---- Parameterized DGP Generators -------------------------------------------

generate_dgp_lognormal <- function(c_shift) {
  p_hurdle_0   <- pnorm(c_shift + 0.5 * X[, 1] - 0.3 * X[, 3])
  p_hurdle_1   <- pnorm(c_shift + 0.5 * X[, 1] - 0.3 * X[, 3] + 0.4 + 0.2 * X[, 1])
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
  
  list(y = Y, true_cate = true_cate, true_ate = mean(true_cate), zero_prop = mean(Y == 0))
}

generate_dgp_gaussian <- function(c_shift) {
  p_hurdle_0   <- pnorm(c_shift + 0.5 * X[, 1] - 0.3 * X[, 3])
  p_hurdle_1   <- pnorm(c_shift + 0.5 * X[, 1] - 0.3 * X[, 3] + 0.4 + 0.2 * X[, 1])
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
  
  list(y = Y, true_cate = true_cate, true_ate = mean(true_cate), zero_prop = mean(Y == 0))
}

generate_dgp_tweedie <- function(c_shift) {
  log_mu0 <- 1.2 + c_shift + 0.8 * X[, 1] - 0.4 * X[, 3]
  log_mu1 <- 1.2 + c_shift + 0.8 * X[, 1] - 0.4 * X[, 3] + 0.6 + 0.3 * X[, 1]
  
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
  
  list(y = Y, true_cate = true_cate, true_ate = mean(true_cate), zero_prop = mean(Y == 0))
}

## ---- Define Scenario Mapping ------------------------------------------------

grid <- data.frame(
  id = 1:15,
  DGP = c(rep("DGP A: Log-Normal Hurdle", 5),
          rep("DGP B: Gaussian Hurdle", 5),
          rep("DGP C: Tweedie Compound", 5)),
  c_shift = c(
    -1.5, -0.5, 0.2, 1.0, 1.8,  # DGP A
    -1.5, -0.5, 0.2, 1.0, 1.8,  # DGP B
    -3.5, -2.0, -0.8, 0.0, 1.0   # DGP C
  ),
  stringsAsFactors = FALSE
)

# Fetch our scenario details
row <- grid[scenario_id, ]
dgp_type <- row$DGP
c_shift <- row$c_shift

# Set scenario isolated seed
set.seed(4200 + scenario_id)

# Generate DGP
if (dgp_type == "DGP A: Log-Normal Hurdle") {
  dgp <- generate_dgp_lognormal(c_shift)
} else if (dgp_type == "DGP B: Gaussian Hurdle") {
  dgp <- generate_dgp_gaussian(c_shift)
} else {
  dgp <- generate_dgp_tweedie(c_shift)
}

cat(sprintf("\n--- [Isolated Scenario %d/15] Starting: %s, c_shift = %.1f, Zero Proportion = %.1f%% ---\n",
            scenario_id, dgp_type, c_shift, 100 * dgp$zero_prop))

# ---- Load Existing Results for baseline and other paths ---------------------
master_res_file <- file.path(RESULTS_DIR, "sensitivity_analysis_results.csv")
if (!file.exists(master_res_file)) {
  stop("Fatal Error: Missing master results file sensitivity_analysis_results.csv.")
}
master_df <- read.csv(master_res_file, stringsAsFactors = FALSE)
scenario_master <- master_df[scenario_id, ]

rmse_linear     <- scenario_master$Linear_RMSE
bias_linear     <- scenario_master$Linear_Bias
coverage_linear <- scenario_master$Linear_Coverage
cor_linear      <- scenario_master$Linear_Correlation

rmse_log     <- scenario_master$Log_RMSE
bias_log     <- scenario_master$Log_Bias
coverage_log <- scenario_master$Log_Coverage
cor_log      <- scenario_master$Log_Correlation

rmse_pathA     <- scenario_master$PathA_RMSE
bias_pathA     <- scenario_master$PathA_Bias
coverage_pathA <- scenario_master$PathA_Coverage
cor_pathA      <- scenario_master$PathA_Correlation

rmse_pathB     <- scenario_master$PathB_RMSE
bias_pathB     <- scenario_master$PathB_Bias
coverage_pathB <- scenario_master$PathB_Coverage
cor_pathB      <- scenario_master$PathB_Correlation

rmse_pathD     <- scenario_master$PathD_RMSE
bias_pathD     <- scenario_master$PathD_Bias
coverage_pathD <- scenario_master$PathD_Coverage
cor_pathD      <- scenario_master$PathD_Correlation

# ----------------------------------------------------
# Model 5: Joint Copula-BCF Path C (Selection model)
# ----------------------------------------------------
cat("  Fitting Joint Copula-BCF (Path C)...\n")
fit_pathC <- pathc_bcf(
    y          = dgp$y,
    z          = Z,
    x_control  = X,
    pihat_sel  = pi_x,
    pihat_out  = NULL,
    nburn      = NBURN,
    nsim       = NSIM,
    nthin      = NTHIN,
    update_interval = 9999
)
cate_draws_pathC <- matrix(0, nrow = NSIM, ncol = N)
for (s in 1:NSIM) {
  eta_b0 <- fit_pathC$sel_con_post[s, ]
  eta_b1 <- fit_pathC$sel_con_post[s, ] + fit_pathC$sel_tau_post[s, ]
  eta_c0 <- fit_pathC$out_con_post[s, ]
  eta_c1 <- fit_pathC$out_con_post[s, ] + fit_pathC$out_tau_post[s, ]
  sig2 <- fit_pathC$sigma_post[s]^2
  bet <- fit_pathC$beta_post[s]
  
  mu0_draw <- exp(eta_c0 + 0.5 * sig2) * pnorm(eta_b0 + bet)
  mu1_draw <- exp(eta_c1 + 0.5 * sig2) * pnorm(eta_b1 + bet)
  cate_draws_pathC[s, ] <- mu1_draw - mu0_draw
}
ate_draws_pathC <- rowMeans(cate_draws_pathC)
cate_est_pathC  <- colMeans(cate_draws_pathC)
cate_ci_pathC   <- apply(cate_draws_pathC, 2, quantile, probs = c(0.025, 0.975))

rmse_pathC     <- sqrt(mean((cate_est_pathC - dgp$true_cate)^2))
bias_pathC     <- mean(cate_est_pathC - dgp$true_cate)
coverage_pathC <- mean(dgp$true_cate >= cate_ci_pathC[1, ] & dgp$true_cate <= cate_ci_pathC[2, ])
cor_pathC      <- cor(cate_est_pathC, dgp$true_cate)

# ----------------------------------------------------
# Model 7: ZIC-BCF-Smear (Best_Path_Gemini)
# ----------------------------------------------------
cat("  Fitting ZIC-BCF-Smear (Best_Path_Gemini)...\n")
fit_smear <- zicbcf_smear(
    y             = dgp$y,
    z             = Z,
    x_control     = X,
    pihat         = pi_x,
    pihat_active  = NULL,
    nburn         = NBURN,
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

cat(sprintf("[Scenario %d] Completed successfully!\n", scenario_id))

# Save scenario checkpoint data frame
out_df <- data.frame(
  DGP = dgp_type,
  c_shift = c_shift,
  Zero_Proportion = dgp$zero_prop,
  True_ATE = dgp$true_ate,
  
  # BCF_Linear
  Linear_RMSE = rmse_linear,
  Linear_Bias = abs(bias_linear),
  Linear_Coverage = coverage_linear,
  Linear_Correlation = cor_linear,
  
  # BCF_Log
  Log_RMSE = rmse_log,
  Log_Bias = abs(bias_log),
  Log_Coverage = coverage_log,
  Log_Correlation = cor_log,
  
  # ZIC-BCF (Path A)
  PathA_RMSE = rmse_pathA,
  PathA_Bias = abs(bias_pathA),
  PathA_Coverage = coverage_pathA,
  PathA_Correlation = cor_pathA,
  
  # Tweedie (Path B)
  PathB_RMSE = rmse_pathB,
  PathB_Bias = abs(bias_pathB),
  PathB_Coverage = coverage_pathB,
  PathB_Correlation = cor_pathB,
  
  # Joint Copula (Path C)
  PathC_RMSE = rmse_pathC,
  PathC_Bias = abs(bias_pathC),
  PathC_Coverage = coverage_pathC,
  PathC_Correlation = cor_pathC,
  
  # Gamma Hurdle (Path D)
  PathD_RMSE = rmse_pathD,
  PathD_Bias = abs(bias_pathD),
  PathD_Coverage = coverage_pathD,
  PathD_Correlation = cor_pathD,
  
  # ZIC-BCF-Smear (Best_Path_Gemini)
  Gemini_RMSE = rmse_smear,
  Gemini_Bias = abs(bias_smear),
  Gemini_Coverage = coverage_smear,
  Gemini_Correlation = cor_smear,
  
  stringsAsFactors = FALSE
)

csv_path <- file.path(RESULTS_DIR, sprintf("scenario_results_%d.csv", scenario_id))
write.csv(out_df, csv_path, row.names = FALSE)
cat(sprintf("[SUCCESS] Scenario %d results saved to: %s\n", scenario_id, csv_path))
