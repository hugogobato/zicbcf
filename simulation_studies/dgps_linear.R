################################################################################
##  Link-linear data-generating processes, extracted verbatim from
##  run_updated_simulation_study.R so that they can be sourced without running
##  the full 100-seed simulation grid.
##
##  Each generator returns the observed data together with the closed-form true
##  potential-outcome means, so any diagnostic that needs a known truth can be
##  run against them.
################################################################################

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

