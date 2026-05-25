# Verification Simulation for Path C Joint Copula-BCF
# Simulates a Heckman selection DGP with selection on unobservables,
# fits the model, performs out-of-sample prediction, and validates recovery.

library(countbcf, lib.loc = "local_lib")

set.seed(42)

cat("\n============================================================\n")
cat("Starting Joint Copula-BCF (Path C) Verification Simulation\n")
cat("============================================================\n\n")

## ---- Scenario Parameters ---------------------------------------------------
N <- 600
P <- 5
rho_true <- 0.5
sigma_true <- 0.5
beta_true <- rho_true * sigma_true  # 0.25
sigma02_true <- sigma_true^2 * (1 - rho_true^2)  # 0.1875
sigma0_true <- sqrt(sigma02_true)  # 0.4330127

cat("True Covariance Parameters:\n")
cat(sprintf("  rho      = %.4f\n", rho_true))
cat(sprintf("  beta     = %.4f\n", beta_true))
cat(sprintf("  sigma    = %.4f\n", sigma_true))
cat(sprintf("  sigma0   = %.4f\n\n", sigma0_true))

## ---- Generate Covariates & Propensity ---------------------------------------
X <- matrix(rnorm(N * P), N, P)
colnames(X) <- paste0("X", 1:P)

# Propensity Score & Treatment Assignment
pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
z <- rbinom(N, 1, pi_x)

## ---- Generate Bivariate Normal Errors ---------------------------------------
# delta_i ~ N(0, 1) (selection error)
# epsilon_i = beta * delta_i + eta_i, eta_i ~ N(0, sigma02) (outcome error)
delta <- rnorm(N, 0, 1)
eta <- rnorm(N, 0, sigma0_true)
epsilon <- beta_true * delta + eta

# Check empirical correlation
emp_corr <- cor(delta, epsilon)
cat(sprintf("Realized empirical error correlation: %.4f\n", emp_corr))

## ---- Hurdle & Intensity Components -----------------------------------------
# Selection (hurdle) prognostic and moderating parts
mu_b_true <- 0.2 + 0.5 * X[, 1] - 0.3 * X[, 3]
tau_b_true <- 0.4 + 0.2 * X[, 1]
eta_b_true <- mu_b_true + z * tau_b_true
W_star_true <- eta_b_true + delta

I <- as.numeric(W_star_true > 0)
active <- I == 1

# Continuous Intensity prognostic and moderating parts
mu_c_true <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4]
tau_c_true <- 0.5 - 0.3 * X[, 2]
eta_c_true <- mu_c_true + z * tau_c_true
V_true <- eta_c_true + epsilon

# Observed semicontinuous outcome
y <- I * exp(V_true)

cat(sprintf("Realized active fraction (hurdle participation): %.2f%%\n", mean(active) * 100))
cat(sprintf("Realized zeros: %.2f%%\n\n", mean(y == 0) * 100))

## ---- Train/Test Split -------------------------------------------------------
# Split 80% train, 20% test
train_idx <- 1:480
test_idx <- 481:600

X_train <- X[train_idx, ]
z_train <- z[train_idx]
y_train <- y[train_idx]
active_train <- active[train_idx]

X_test <- X[test_idx, ]
z_test <- z[test_idx]
y_test <- y[test_idx]
active_test <- active[test_idx]

# True response potential outcomes (training scale)
# E[Y(z) | X] = exp(eta_c(z) + sigma^2 / 2) * Phi(eta_b(z) + beta)
calc_true_po <- function(X_mat, z_val) {
  mu_b <- 0.2 + 0.5 * X_mat[, 1] - 0.3 * X_mat[, 3]
  tau_b <- 0.4 + 0.2 * X_mat[, 1]
  eta_b <- mu_b + z_val * tau_b
  
  mu_c <- 1.5 + 0.8 * X_mat[, 2] + 0.4 * X_mat[, 4]
  tau_c <- 0.5 - 0.3 * X_mat[, 2]
  eta_c <- mu_c + z_val * tau_c
  
  exp(eta_c + 0.5 * sigma_true^2) * pnorm(eta_b + beta_true)
}

true_mu0_train <- calc_true_po(X_train, 0)
true_mu1_train <- calc_true_po(X_train, 1)
true_cate_train <- true_mu1_train - true_mu0_train

true_mu0_test <- calc_true_po(X_test, 0)
true_mu1_test <- calc_true_po(X_test, 1)
true_cate_test <- true_mu1_test - true_mu0_test

## ---- Fit Model with Out-Of-Sample Prediction --------------------------------
cat("Fitting Joint Copula-BCF with out-of-sample predictions...\n")
t0 <- Sys.time()
fit <- pathc_bcf(
  y = y_train, z = z_train,
  x_control = X_train,
  nburn = 400, nsim = 600, nthin = 1,
  update_interval = 200,
  x_control_est = X_test,
  z_est = z_test
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Fitting finished in %.1f seconds.\n\n", elapsed))

## =================== COVARIATION PARAMETERS RECOVERY ===================
cat("=================== COVARIANCE PARAMETERS ===================\n")
get_stats <- function(draws, true_val) {
  p_mean <- mean(draws)
  p_sd   <- sd(draws)
  ci     <- quantile(draws, c(0.025, 0.975))
  sprintf("True = %6.3f | Post Mean = %6.3f | Post SD = %6.3f | 95%% CI = [%6.3f, %6.3f]",
          true_val, p_mean, p_sd, ci[1], ci[2])
}

cat("rho   : ", get_stats(fit$rho_post, rho_true), "\n")
cat("beta  : ", get_stats(fit$beta_post, beta_true), "\n")
cat("sigma : ", get_stats(fit$sigma_post, sigma_true), "\n")
cat("sigma0: ", get_stats(fit$sigma0_post, sigma0_true), "\n\n")

## =================== LATENT VARIABLES RECOVERY ===================
cat("=================== LATENT VARIABLES RECOVERY ===================\n")
# In-sample latent variables
W_post_mean <- colMeans(fit$W_post)
V_post_mean <- colMeans(fit$V_post)

corr_W <- cor(W_post_mean, W_star_true[train_idx])
# For V, we only evaluate correlation on censored (zero) units, since active units are fixed at observed y
censored_idx <- which(y_train == 0)
corr_V <- cor(V_post_mean[censored_idx], V_true[train_idx][censored_idx])

cat(sprintf("Correlation between posterior mean W* and true latent W* (all): %.4f\n", corr_W))
cat(sprintf("Correlation between posterior mean V and true latent V (censored only): %.4f\n\n", corr_V))

## =================== IN-SAMPLE ACCURACY ===================
cat("=================== IN-SAMPLE RESPONSE ACCURACY ===================\n")
# Estimate potential outcomes in-sample
# For each draw, calculate the response scale potential outcomes
nsim <- nrow(fit$out_con_post)
n_train <- ncol(fit$out_con_post)

mu0_post_train <- matrix(0, nsim, n_train)
mu1_post_train <- matrix(0, nsim, n_train)

for (s in 1:nsim) {
  # active log-outcome SD and mean
  muy <- fit$muy_active
  sdy <- fit$sdy_active
  
  # Selection probability parameter
  eta_b0 <- fit$sel_con_post[s, ]
  eta_b1 <- fit$sel_con_post[s, ] + fit$sel_mod_post[s, ]
  
  # Log intensity parameter (un-scaled log-normal parameters)
  eta_c0 <- fit$out_con_post[s, ]
  eta_c1 <- fit$out_con_post[s, ] + fit$out_mod_post[s, ]
  
  # Total outcome variance (original scale)
  sig2 <- fit$sigma_post[s]^2
  
  # beta regression coefficient (original scale)
  bet <- fit$beta_post[s]
  
  # Expected values on original physical scale
  mu0_post_train[s, ] <- exp(eta_c0 + 0.5 * sig2) * pnorm(eta_b0 + bet)
  mu1_post_train[s, ] <- exp(eta_c1 + 0.5 * sig2) * pnorm(eta_b1 + bet)
}

cate_post_train <- mu1_post_train - mu0_post_train
cate_est_train <- colMeans(cate_post_train)
rmse_cate_train <- sqrt(mean((cate_est_train - true_cate_train)^2))
bias_cate_train <- mean(cate_est_train - true_cate_train)

cat(sprintf("In-sample CATE RMSE : %.4f\n", rmse_cate_train))
cat(sprintf("In-sample CATE Bias : %+.4f\n", bias_cate_train))
cat(sprintf("In-sample True ATE  : %.4f | Est ATE = %.4f (95%% CI = [%.4f, %.4f])\n\n",
            mean(true_cate_train), mean(cate_post_train),
            quantile(rowMeans(cate_post_train), 0.025),
            quantile(rowMeans(cate_post_train), 0.975)))

## =================== OUT-OF-SAMPLE ACCURACY ===================
cat("=================== OUT-OF-SAMPLE RESPONSE ACCURACY ===================\n")
# Estimate potential outcomes out-of-sample
n_test <- ncol(fit$out_con_est_post)

mu0_post_test <- matrix(0, nsim, n_test)
mu1_post_test <- matrix(0, nsim, n_test)

for (s in 1:nsim) {
  # Selection probability parameter (out of sample draws)
  eta_b0 <- fit$sel_con_est_post[s, ]
  eta_b1 <- fit$sel_con_est_post[s, ] + fit$sel_mod_est_post[s, ]
  
  # Log intensity parameter (un-scaled log-normal parameters)
  eta_c0 <- fit$out_con_est_post[s, ]
  eta_c1 <- fit$out_con_est_post[s, ] + fit$out_mod_est_post[s, ]
  
  # Total outcome variance (original scale)
  sig2 <- fit$sigma_post[s]^2
  
  # beta regression coefficient (original scale)
  bet <- fit$beta_post[s]
  
  # Expected values on original physical scale
  mu0_post_test[s, ] <- exp(eta_c0 + 0.5 * sig2) * pnorm(eta_b0 + bet)
  mu1_post_test[s, ] <- exp(eta_c1 + 0.5 * sig2) * pnorm(eta_b1 + bet)
}

cate_post_test <- mu1_post_test - mu0_post_test
cate_est_test <- colMeans(cate_post_test)
rmse_cate_test <- sqrt(mean((cate_est_test - true_cate_test)^2))
bias_cate_test <- mean(cate_est_test - true_cate_test)

cat(sprintf("Out-of-sample CATE RMSE : %.4f\n", rmse_cate_test))
cat(sprintf("Out-of-sample CATE Bias : %+.4f\n", bias_cate_test))
cat(sprintf("Out-of-sample True ATE  : %.4f | Est ATE = %.4f (95%% CI = [%.4f, %.4f])\n",
            mean(true_cate_test), mean(cate_post_test),
            quantile(rowMeans(cate_post_test), 0.025),
            quantile(rowMeans(cate_post_test), 0.975)))
cat("============================================================\n")
