################################################################################
##  Semicontinuous BCF Extended Study: Single-Seed Process Executor
##
##  Fits BCF-Linear, ZIC-BCF (Path A), and ZIC-BCF-Smear on a single seed
##  for a given DGP to ensure clean memory release upon completion.
##
##  Usage:
##    Rscript simulation_studies/reduced_simulations/run_5seed_single.R <DGP> <seed>
################################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Error: Must specify DGP ('A', 'B', or 'C') and seed (integer) as arguments.")
}

dgp_type <- args[1]
seed_num <- as.integer(args[2])

if (!(dgp_type %in% c("A", "B", "C"))) {
  stop("Error: DGP must be 'A', 'B', or 'C'")
}
if (is.na(seed_num)) {
  stop("Error: Seed must be a valid integer")
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

# Generate dataset-specific covariates and assignments under seed
set.seed(seed_num * 1000 + 42)
X <- matrix(rnorm(N * P), N, P)
colnames(X) <- paste0("X", 1:P)

# Propensity score (Confounded)
pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
Z    <- rbinom(N, 1, pi_x)

## ---- Parameterized DGP Generators -------------------------------------------

if (dgp_type == "A") {
  # DGP A: Log-Normal Hurdle
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
  
  true_hurdle_cate <- p_hurdle_1 - p_hurdle_0
  dgp_name <- "DGP A: Log-Normal Hurdle"
  
} else if (dgp_type == "B") {
  # DGP B: Gaussian Hurdle
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
  
  true_hurdle_cate <- p_hurdle_1 - p_hurdle_0
  dgp_name <- "DGP B: Gaussian Hurdle"
  
} else {
  # DGP C: Tweedie Compound Poisson-Gamma (p=1.5)
  log_mu0 <- 1.2 + 0.8 * X[, 1] - 0.4 * X[, 3]
  log_mu1 <- 1.2 + 0.8 * X[, 1] - 0.4 * X[, 3] + 0.6 + 0.3 * X[, 1]
  
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
  
  p0_hurdle_true <- 1 - exp(-lambda0_true)
  p1_hurdle_true <- 1 - exp(-lambda1_true)
  true_hurdle_cate <- p1_hurdle_true - p0_hurdle_true
  dgp_name <- "DGP C: Tweedie Compound"
}

true_ate <- mean(true_cate)
true_hurdle_ate <- mean(true_hurdle_cate)

cat(sprintf("\n--- Fitting models for %s | Seed %d ---\n", dgp_name, seed_num))

# Helper to pack standard CATE metrics
calc_cate_metrics <- function(cate_draws, true_c, ate_draws) {
  cate_est <- colMeans(cate_draws)
  cate_ci  <- apply(cate_draws, 2, quantile, probs = c(0.025, 0.975))
  
  rmse <- sqrt(mean((cate_est - true_c)^2))
  bias <- mean(cate_est - true_c)
  coverage <- mean(true_c >= cate_ci[1, ] & true_c <= cate_ci[2, ])
  correlation <- cor(cate_est, true_c)
  ci_length <- mean(cate_ci[2, ] - cate_ci[1, ])
  est_ate_mean <- mean(ate_draws)
  
  list(rmse=rmse, bias=bias, coverage=coverage, correlation=correlation, ci_length=ci_length, est_ate_mean=est_ate_mean)
}

# ----------------------------------------------------
# Model 1: BCF-Linear
# ----------------------------------------------------
cat("  1. Fitting BCF-Linear...\n")
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
cat("  2. Fitting ZIC-BCF (Path A)...\n")
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

# Hurdle metrics for Path A
p0_draws_pathA <- pnorm(fit_pathA$mu_b)
p1_draws_pathA <- pnorm(fit_pathA$mu_b + fit_pathA$tau_b)
hurdle_cate_draws_pathA <- p1_draws_pathA - p0_draws_pathA
hurdle_ate_draws_pathA  <- rowMeans(hurdle_cate_draws_pathA)
m_hurdle_pathA <- calc_cate_metrics(hurdle_cate_draws_pathA, true_hurdle_cate, hurdle_ate_draws_pathA)

# ----------------------------------------------------
# Model 3: ZIC-BCF-Smear (Best_Path_Gemini)
# ----------------------------------------------------
cat("  3. Fitting ZIC-BCF-Smear (Best_Path_Gemini)...\n")
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

# Hurdle metrics for Smear (same as Path A probit stage draws)
p0_draws_smear <- pnorm(fit_smear$mu_b)
p1_draws_smear <- pnorm(fit_smear$mu_b + fit_smear$tau_b)
hurdle_cate_draws_smear <- p1_draws_smear - p0_draws_smear
hurdle_ate_draws_smear  <- rowMeans(hurdle_cate_draws_smear)
m_hurdle_smear <- calc_cate_metrics(hurdle_cate_draws_smear, true_hurdle_cate, hurdle_ate_draws_smear)

# ----------------------------------------------------
# Save Results CSV
# ----------------------------------------------------
df_results <- data.frame(
  DGP = dgp_name,
  Seed = seed_num,
  True_ATE = true_ate,
  True_Hurdle_ATE = true_hurdle_ate,
  
  # Standard CATE Metrics
  Linear_CATE_RMSE        = m_linear$rmse,
  Linear_CATE_Abs_Bias    = abs(m_linear$bias),
  Linear_CATE_Coverage    = m_linear$coverage,
  Linear_CATE_Correlation = m_linear$correlation,
  Linear_CATE_CI_Length   = m_linear$ci_length,
  Linear_Est_ATE          = m_linear$est_ate_mean,
  
  PathA_CATE_RMSE        = m_pathA$rmse,
  PathA_CATE_Abs_Bias    = abs(m_pathA$bias),
  PathA_CATE_Coverage    = m_pathA$coverage,
  PathA_CATE_Correlation = m_pathA$correlation,
  PathA_CATE_CI_Length   = m_pathA$ci_length,
  PathA_Est_ATE          = m_pathA$est_ate_mean,
  
  Smear_CATE_RMSE        = m_smear$rmse,
  Smear_CATE_Abs_Bias    = abs(m_smear$bias),
  Smear_CATE_Coverage    = m_smear$coverage,
  Smear_CATE_Correlation = m_smear$correlation,
  Smear_CATE_CI_Length   = m_smear$ci_length,
  Smear_Est_ATE          = m_smear$est_ate_mean,
  
  # Hurdle CATE Metrics (Linear is NA)
  Linear_Hurdle_RMSE        = NA,
  Linear_Hurdle_Abs_Bias    = NA,
  Linear_Hurdle_Coverage    = NA,
  Linear_Hurdle_Correlation = NA,
  Linear_Hurdle_CI_Length   = NA,
  Linear_Est_Hurdle_ATE     = NA,
  
  PathA_Hurdle_RMSE        = m_hurdle_pathA$rmse,
  PathA_Hurdle_Abs_Bias    = abs(m_hurdle_pathA$bias),
  PathA_Hurdle_Coverage    = m_hurdle_pathA$coverage,
  PathA_Hurdle_Correlation = m_hurdle_pathA$correlation,
  PathA_Hurdle_CI_Length   = m_hurdle_pathA$ci_length,
  PathA_Est_Hurdle_ATE     = m_hurdle_pathA$est_ate_mean,
  
  Smear_Hurdle_RMSE        = m_hurdle_smear$rmse,
  Smear_Hurdle_Abs_Bias    = abs(m_hurdle_smear$bias),
  Smear_Hurdle_Coverage    = m_hurdle_smear$coverage,
  Smear_Hurdle_Correlation = m_hurdle_smear$correlation,
  Smear_Hurdle_CI_Length   = m_hurdle_smear$ci_length,
  Smear_Est_Hurdle_ATE     = m_hurdle_smear$est_ate_mean,
  
  stringsAsFactors = FALSE
)

csv_filename <- sprintf("checkpoint_5seed_DGP_%s_seed_%d.csv", dgp_type, seed_num)
csv_filepath <- file.path(RESULTS_DIR, csv_filename)
write.csv(df_results, csv_filepath, row.names = FALSE)

cat(sprintf("[SUCCESS] Completed DGP %s Seed %d. Saved to %s\n", dgp_type, seed_num, csv_filepath))
