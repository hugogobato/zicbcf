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

# Propensity Score & Treatment Assignment (Confounded)
pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
Z    <- rbinom(N, 1, pi_x)

## ---- DGP generators ---------------------------------------------------------
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
  
  list(y = Y, true_cate = true_cate, true_ate = mean(true_cate), name = "Log-Normal Hurdle (DGP A)")
}

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
  
  list(y = Y, true_cate = true_cate, true_ate = mean(true_cate), name = "Gaussian Hurdle (DGP B)")
}

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
  
  list(y = Y, true_cate = true_cate, true_ate = mean(true_cate), name = "Tweedie Compound (DGP C)")
}

## ---- Evaluator for Path D only ---------------------------------------------
eval_pathd <- function(dgp, nburn) {
  cat(sprintf("\n=== Evaluating %s with nburn=%d ===\n", dgp$name, nburn))
  
  # Step 1: Fit hurdle stage using zicbcf_pathA to get exact same probit hurdle draws
  cat("  Fitting Hurdle Probit BCF stage...\n")
  fit_hurdle <- zicbcf_pathA(
      y             = dgp$y,
      z             = Z,
      x_control     = X,
      pihat         = pi_x,
      pihat_active  = NULL,
      nburn         = nburn,
      nsim          = NSIM,
      nthin         = NTHIN,
      update_interval = 9999
  )
  mu_b_all  <- fit_hurdle$mu_b
  tau_b_all <- fit_hurdle$tau_b
  
  p0_hurdle <- pnorm(mu_b_all)
  p1_hurdle <- pnorm(mu_b_all + tau_b_all)
  
  # Step 2: Fit updated Gamma intensity stage
  cat("  Fitting Intensity Log-Linear Gamma BCF stage...\n")
  fit_intensity <- pathd_gammabcf(
      y             = dgp$y,
      z             = Z,
      x_control     = X,
      pihat_pos     = NULL, # SPA automatically estimated
      nburn         = nburn,
      nsim          = NSIM,
      thin          = NTHIN,
      update_interval = 9999,
      return_trees  = TRUE
  )
  
  fit_intensity$control_fit$scale <- 1.0
  fit_intensity$control_fit$shift <- 0.0
  fit_intensity$moderate_fit$scale <- 1.0
  fit_intensity$moderate_fit$shift <- 0.0
  
  # Correct Prediction using SPA propensity score column
  X_c_all <- cbind(X, pi_x)
  mu_f_all  <- -get_forest_fit(fit_intensity$control_fit, X_c_all)
  tau_f_all <- -get_forest_fit(fit_intensity$moderate_fit, X)
  
  lambda_0_all <- exp(mu_f_all)
  lambda_1_all <- exp(mu_f_all + tau_f_all)
  
  cate_draws_pathD <- matrix(0, nrow = NSIM, ncol = N)
  for (s in 1:NSIM) {
    mu0_draw <- p0_hurdle[s, ] * lambda_0_all[s, ]
    mu1_draw <- p1_hurdle[s, ] * lambda_1_all[s, ]
    cate_draws_pathD[s, ] <- mu1_draw - mu0_draw
  }
  ate_draws_pathD <- rowMeans(cate_draws_pathD)
  
  cate_est_pathD  <- colMeans(cate_draws_pathD)
  cate_ci_pathD   <- apply(cate_draws_pathD, 2, quantile, probs = c(0.025, 0.975))
  
  rmse_pathD     <- sqrt(mean((cate_est_pathD - dgp$true_cate)^2))
  bias_pathD     <- mean(cate_est_pathD - dgp$true_cate)
  coverage_pathD <- mean(dgp$true_cate >= cate_ci_pathD[1, ] & dgp$true_cate <= cate_ci_pathD[2, ])
  cor_pathD      <- cor(cate_est_pathD, dgp$true_cate)
  
  cat(sprintf("Estimated shape parameter kappa_c (mean posterior draw): %.4f\n", mean(fit_intensity$kappa)))
  cat(sprintf("Path D updated results: Est ATE Mean = %.4f, Est ATE SD = %.4f, CATE RMSE = %.4f, CATE Abs Bias = %.4f, CATE 95%% Coverage = %.2f%%, CATE Correlation = %.4f\n",
              mean(ate_draws_pathD), sd(ate_draws_pathD), rmse_pathD, abs(bias_pathD), 100 * coverage_pathD, cor_pathD))
  
  return(c(mean(ate_draws_pathD), sd(ate_draws_pathD), rmse_pathD, abs(bias_pathD), coverage_pathD, cor_pathD))
}

# Run the evaluations
dgp_a <- generate_dgp_lognormal()
res_a <- eval_pathd(dgp_a, nburn = 500)

dgp_b <- generate_dgp_gaussian()
res_b <- eval_pathd(dgp_b, nburn = 500)

dgp_c <- generate_dgp_tweedie()
res_c_500 <- eval_pathd(dgp_c, nburn = 500)
res_c_2000 <- eval_pathd(dgp_c, nburn = 2000)

cat("\n=================== FINAL SUMMARY TABLE ===================\n")
df_summary <- data.frame(
  DGP = c("DGP A: Log-Normal Hurdle", "DGP B: Gaussian Hurdle", "DGP C: Tweedie Compound (nburn=500)", "DGP C: Tweedie Compound (nburn=2000)"),
  Est_ATE_Mean = c(res_a[1], res_b[1], res_c_500[1], res_c_2000[1]),
  Est_ATE_SD = c(res_a[2], res_b[2], res_c_500[2], res_c_2000[2]),
  CATE_RMSE = c(res_a[3], res_b[3], res_c_500[3], res_c_2000[3]),
  CATE_Abs_Bias = c(res_a[4], res_b[4], res_c_500[4], res_c_2000[4]),
  CATE_95_Coverage = c(res_a[5], res_b[5], res_c_500[5], res_c_2000[5]),
  CATE_Correlation = c(res_a[6], res_b[6], res_c_500[6], res_c_2000[6])
)
print(df_summary, digits = 4)
