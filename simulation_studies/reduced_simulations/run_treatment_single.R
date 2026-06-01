################################################################################
##  Semicontinuous BCF Treatment Sensitivity: Single Scenario Executor
##
##  Fits BCF-Linear, ZIC-BCF (Path A), and ZIC-BCF-Smear for a single scenario
##  (1 to 15) representing a combination of DGP and treatment magnitude multiplier k.
##
##  Usage:
##    Rscript simulation_studies/reduced_simulations/run_treatment_single.R <scenario_id>
################################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Error: Must specify scenario ID (1 to 15) as an argument.")
}
scenario_id <- as.integer(args[1])
if (is.na(scenario_id) || scenario_id < 1 || scenario_id > 15) {
  stop("Error: Scenario ID must be an integer between 1 and 15.")
}

library(countbcf, lib.loc = "local_lib")

# Output directory setup
RESULTS_DIR <- "simulation_studies/reduced_simulations/results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

## ---- Global Simulation Parameters --------------------------------------------
N     <- 1000
P     <- 5
NBURN <- 2000  # High burn-in to ensure perfect convergence
NSIM  <- 1000
NTHIN <- 1

# Scenario grid definition
grid <- data.frame(
  id = 1:15,
  DGP = c(rep("DGP A: Log-Normal Hurdle", 5),
          rep("DGP B: Gaussian Hurdle", 5),
          rep("DGP C: Tweedie Compound", 5)),
  dgp_code = c(rep("A", 5), rep("B", 5), rep("C", 5)),
  k = rep(c(0.0, 0.5, 1.0, 1.5, 2.0), 3),
  stringsAsFactors = FALSE
)

scenario <- grid[scenario_id, ]
dgp_code <- scenario$dgp_code
dgp_name <- scenario$DGP
k        <- scenario$k

# Set scenario-specific seed
set.seed(8800 + scenario_id)

# Generate covariates and confounded assignments
X <- matrix(rnorm(N * P), N, P)
colnames(X) <- paste0("X", 1:P)
pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
Z    <- rbinom(N, 1, pi_x)

## ---- Parameterized DGP Generators with Treatment Multiplier k ----------------

if (dgp_code == "A") {
  # DGP A: Log-Normal Hurdle (target ~60% zero-inflation)
  p_hurdle_0   <- pnorm(-0.3 + 0.5 * X[, 1] - 0.3 * X[, 3])
  p_hurdle_1   <- pnorm(-0.3 + 0.5 * X[, 1] - 0.3 * X[, 3] + k * (0.4 + 0.2 * X[, 1]))
  p_hurdle_obs <- ifelse(Z == 1, p_hurdle_1, p_hurdle_0)
  I <- rbinom(N, 1, p_hurdle_obs)
  
  mu_c_0     <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4]
  mu_c_1     <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4] + k * (0.5 - 0.3 * X[, 2])
  sigma_true <- 0.5
  
  y_pos_0   <- exp(mu_c_0 + rnorm(N, 0, sigma_true))
  y_pos_1   <- exp(mu_c_1 + rnorm(N, 0, sigma_true))
  y_pos_obs <- ifelse(Z == 1, y_pos_1, y_pos_0)
  Y <- I * y_pos_obs
  
  true_mu0  <- p_hurdle_0 * exp(mu_c_0 + 0.5 * sigma_true^2)
  true_mu1  <- p_hurdle_1 * exp(mu_c_1 + 0.5 * sigma_true^2)
  true_cate <- true_mu1 - true_mu0
  
} else if (dgp_code == "B") {
  # DGP B: Gaussian Hurdle (target ~60% zero-inflation)
  p_hurdle_0   <- pnorm(-0.3 + 0.5 * X[, 1] - 0.3 * X[, 3])
  p_hurdle_1   <- pnorm(-0.3 + 0.5 * X[, 1] - 0.3 * X[, 3] + k * (0.4 + 0.2 * X[, 1]))
  p_hurdle_obs <- ifelse(Z == 1, p_hurdle_1, p_hurdle_0)
  I <- rbinom(N, 1, p_hurdle_obs)
  
  mu_c_0     <- 6.0 + 1.2 * X[, 2] + 0.6 * X[, 4]
  mu_c_1     <- 6.0 + 1.2 * X[, 2] + 0.6 * X[, 4] + k * (1.5 - 0.5 * X[, 2])
  sigma_true <- 1.0
  
  y_pos_0   <- mu_c_0 + rnorm(N, 0, sigma_true)
  y_pos_1   <- mu_c_1 + rnorm(N, 0, sigma_true)
  y_pos_obs <- ifelse(Z == 1, y_pos_1, y_pos_0)
  Y <- I * y_pos_obs
  
  true_mu0  <- p_hurdle_0 * mu_c_0
  true_mu1  <- p_hurdle_1 * mu_c_1
  true_cate <- true_mu1 - true_mu0
  
} else {
  # DGP C: Tweedie Compound Poisson-Gamma (p=1.5) (target ~60% zero-inflation via -1.0 shift)
  log_mu0 <- 1.2 - 1.0 + 0.8 * X[, 1] - 0.4 * X[, 3]
  log_mu1 <- 1.2 - 1.0 + 0.8 * X[, 1] - 0.4 * X[, 3] + k * (0.6 + 0.3 * X[, 1])
  
  mu0_true  <- exp(log_mu0)
  mu1_true  <- exp(log_mu1)
  true_cate <- mu1_true - mu0_true
  
  mu_true  <- ifelse(Z == 1, mu1_true, mu0_true)
  phi_true <- 1.5
  
  lambda0_true <- 2 * sqrt(mu0_true) / phi_true
  lambda1_true <- 2 * sqrt(mu1_true) / phi_true
  
  lambda_true <- ifelse(Z == 1, lambda1_true, lambda0_true)
  N_latent    <- rpois(N, lambda_true)
  
  gamma_true0 <- 0.5 * phi_true * sqrt(mu0_true)
  gamma_true1 <- 0.5 * phi_true * sqrt(mu1_true)
  gamma_true  <- ifelse(Z == 1, gamma_true1, gamma_true0)
  
  Y <- rep(0, N)
  for (i in 1:N) {
    if (N_latent[i] > 0) {
      Y[i] <- rgamma(1, shape = N_latent[i], scale = gamma_true[i])
    }
  }
}

true_ate <- mean(true_cate)
zero_prop <- mean(Y == 0)

cat(sprintf("\n--- [Scenario %d/15] Starting: %s, k = %.1f, True ATE = %.3f, Zeros = %.1f%% ---\n",
            scenario_id, dgp_name, k, true_ate, 100 * zero_prop))

# Helper to pack CATE metrics
calc_cate_metrics <- function(cate_draws, true_c, ate_draws) {
  cate_est <- colMeans(cate_draws)
  cate_ci  <- apply(cate_draws, 2, quantile, probs = c(0.025, 0.975))
  
  rmse <- sqrt(mean((cate_est - true_c)^2))
  bias <- mean(cate_est - true_c)
  coverage <- mean(true_c >= cate_ci[1, ] & true_c <= cate_ci[2, ])
  correlation <- cor(cate_est, true_c)
  if (is.na(correlation)) correlation <- 0.0 # handle zero treatment effect case k=0
  
  list(rmse = rmse, bias = bias, coverage = coverage, correlation = correlation, est_ate = mean(ate_draws))
}

# ----------------------------------------------------
# Model 1: BCF-Linear
# ----------------------------------------------------
cat("  Fitting BCF-Linear...\n")
fit_linear <- bcf_continuous_linear(
    y          = Y,
    z          = Z,
    x_control  = X,
    x_moderate = X,
    zhat       = pi_x,
    nburn      = NBURN,
    nsim       = NSIM,
    nthin      = NTHIN,
    update_interval = 9999
)
cate_draws_linear <- get_forest_fit(fit_linear$moderate_fit, X)
ate_draws_linear  <- rowMeans(cate_draws_linear)
m_linear <- calc_cate_metrics(cate_draws_linear, true_cate, ate_draws_linear)

# ----------------------------------------------------
# Model 2: ZIC-BCF (Path A)
# ----------------------------------------------------
cat("  Fitting ZIC-BCF (Path A)...\n")
fit_pathA <- zicbcf_pathA(
    y             = Y,
    z             = Z,
    x_control     = X,
    pihat         = pi_x,
    pihat_active  = NULL,
    nburn         = NBURN,
    nsim          = NSIM,
    nthin         = NTHIN,
    update_interval = 9999
)
m_pathA <- calc_cate_metrics(fit_pathA$cate, true_cate, fit_pathA$ate)

# ----------------------------------------------------
# Model 3: ZIC-BCF-Smear (Best_Path_Gemini)
# ----------------------------------------------------
cat("  Fitting ZIC-BCF-Smear (Best_Path_Gemini)...\n")
fit_smear <- zicbcf_smear(
    y             = Y,
    z             = Z,
    x_control     = X,
    pihat         = pi_x,
    pihat_active  = NULL,
    nburn         = NBURN,
    nsim          = NSIM,
    nthin         = NTHIN,
    update_interval = 9999
)
m_smear <- calc_cate_metrics(fit_smear$cate, true_cate, fit_smear$ate)

# Save checkpoint
out_df <- data.frame(
  Scenario_ID = scenario_id,
  DGP = dgp_name,
  k = k,
  True_ATE = true_ate,
  Zero_Proportion = zero_prop,
  
  # BCF_Linear
  Linear_RMSE = m_linear$rmse,
  Linear_Bias = abs(m_linear$bias),
  Linear_Coverage = m_linear$coverage,
  Linear_Correlation = m_linear$correlation,
  Linear_Est_ATE = m_linear$est_ate,
  
  # ZIC_BCF (Path A)
  PathA_RMSE = m_pathA$rmse,
  PathA_Bias = abs(m_pathA$bias),
  PathA_Coverage = m_pathA$coverage,
  PathA_Correlation = m_pathA$correlation,
  PathA_Est_ATE = m_pathA$est_ate,
  
  # ZIC_BCF_Smear (Best_Path_Gemini)
  Gemini_RMSE = m_smear$rmse,
  Gemini_Bias = abs(m_smear$bias),
  Gemini_Coverage = m_smear$coverage,
  Gemini_Correlation = m_smear$correlation,
  Gemini_Est_ATE = m_smear$est_ate,
  
  stringsAsFactors = FALSE
)

csv_filepath <- file.path(RESULTS_DIR, sprintf("scenario_treatment_results_%d.csv", scenario_id))
write.csv(out_df, csv_filepath, row.names = FALSE)
cat(sprintf("[SUCCESS] Scenario %d results saved to: %s\n", scenario_id, csv_filepath))
