################################################################################
##  Full Semicontinuous BCF Simulation Study
##
##  This script evaluates six methods across three semicontinuous DGPs:
##    - DGP A: Log-Normal Hurdle DGP (proposal standard)
##    - DGP B: Gaussian (Normal) Hurdle DGP (continuous component is normal)
##    - DGP C: Intrinsic Zero-Inflated DGP (Tweedie compound Poisson-Gamma, p=1.5)
##
##  Methods evaluated:
##    1. BCF-Linear: Standard BCF fitted directly on the raw semicontinuous Y
##    2. BCF-Log: Standard BCF fitted on log(Y + 1) and re-transformed
##    3. ZIC-BCF (Path A): Two-part hurdle BCF with SPA propensity adjustment
##    4. Tweedie BCF (Path B): Single-stage compound Poisson-Gamma BCF (with fixed power p=1.5) using exact GIG-conjugacy
##    5. Joint Copula-BCF (Path C): Hurdle selection model with copula correlation
##    6. Gamma Hurdle BCF (Path D): Two-part Gamma BCF with GIG conjugacy and SPA
################################################################################

library(countbcf, lib.loc = "local_lib")

set.seed(42)

## ---- Global specs -----------------------------------------------------------
N     <- 1000
P     <- 5
NBURN <- 500
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
  # True log-mean parameter
  log_mu0 <- 1.2 + 0.8 * X[, 1] - 0.4 * X[, 3]
  log_mu1 <- 1.2 + 0.8 * X[, 1] - 0.4 * X[, 3] + 0.6 + 0.3 * X[, 1]
  
  mu0_true  <- exp(log_mu0)
  mu1_true  <- exp(log_mu1)
  true_cate <- mu1_true - mu0_true
  
  mu_true  <- ifelse(Z == 1, mu1_true, mu0_true)
  phi_true <- 1.5
  
  # Simulation via latent Poisson and continuous Gamma sizes
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

## ---- Fitting & Evaluation Loop ----------------------------------------------
run_eval <- function(dgp) {
  cat(sprintf("\n========================================================================\n"))
  cat(sprintf("=== Evaluating DGP: %s (True ATE = %.4f) ===\n", dgp$name, dgp$true_ate))
  cat(sprintf("========================================================================\n"))
  
  # ----------------------------------------------------
  # Model 1: BCF-Linear (fitting on raw Y)
  # ----------------------------------------------------
  cat("\n  -> Fitting BCF-Linear...\n")
  fit_linear <- bcf_continuous_linear(
      y          = dgp$y,
      z          = Z,
      x_control  = X,
      x_moderate = X,
      zhat       = pi_x,
      nburn      = NBURN,
      nsim       = NSIM,
      nthin      = NTHIN,
      update_interval = 500
  )
  
  cate_draws_linear <- get_forest_fit(fit_linear$moderate_fit, X)
  ate_draws_linear  <- rowMeans(cate_draws_linear)
  
  cate_est_linear   <- colMeans(cate_draws_linear)
  cate_ci_linear    <- apply(cate_draws_linear, 2, quantile, probs = c(0.025, 0.975))
  
  rmse_linear     <- sqrt(mean((cate_est_linear - dgp$true_cate)^2))
  bias_linear     <- mean(cate_est_linear - dgp$true_cate)
  coverage_linear <- mean(dgp$true_cate >= cate_ci_linear[1, ] & dgp$true_cate <= cate_ci_linear[2, ])
  cor_linear      <- cor(cate_est_linear, dgp$true_cate)
  
  # ----------------------------------------------------
  # Model 2: BCF-Log (fitting on log(Y+1))
  # ----------------------------------------------------
  cat("\n  -> Fitting BCF-Log...\n")
  Y_log <- log(dgp$y + 1)
  muy_log <- mean(Y_log)
  fit_log <- bcf_continuous_linear(
      y          = Y_log,
      z          = Z,
      x_control  = X,
      x_moderate = X,
      zhat       = pi_x,
      nburn      = NBURN,
      nsim       = NSIM,
      nthin      = NTHIN,
      update_interval = 500
  )
  
  mu_post_log   <- muy_log + get_forest_fit(fit_log$control_fit, X)
  tau_post_log  <- get_forest_fit(fit_log$moderate_fit, X)
  sigma_post_log <- fit_log$sigma
  
  # Re-transform: Y = exp(W) - 1
  cate_draws_log <- matrix(0, nrow = NSIM, ncol = N)
  for (s in 1:NSIM) {
    mu0_draw <- exp(mu_post_log[s, ] + 0.5 * sigma_post_log[s]^2) - 1
    mu1_draw <- exp(mu_post_log[s, ] + tau_post_log[s, ] + 0.5 * sigma_post_log[s]^2) - 1
    cate_draws_log[s, ] <- mu1_draw - mu0_draw
  }
  ate_draws_log <- rowMeans(cate_draws_log)
  
  cate_est_log  <- colMeans(cate_draws_log)
  cate_ci_log   <- apply(cate_draws_log, 2, quantile, probs = c(0.025, 0.975))
  
  rmse_log     <- sqrt(mean((cate_est_log - dgp$true_cate)^2))
  bias_log     <- mean(cate_est_log - dgp$true_cate)
  coverage_log <- mean(dgp$true_cate >= cate_ci_log[1, ] & dgp$true_cate <= cate_ci_log[2, ])
  cor_log      <- cor(cate_est_log, dgp$true_cate)
  
  # ----------------------------------------------------
  # Model 3: ZIC-BCF Path A (Two-part BCF with SPA)
  # ----------------------------------------------------
  cat("\n  -> Fitting ZIC-BCF (Path A)...\n")
  fit_pathA <- zicbcf_pathA(
      y             = dgp$y,
      z             = Z,
      x_control     = X,
      pihat         = pi_x,
      pihat_active  = NULL, # estimated automatically via SPA
      nburn         = NBURN,
      nsim          = NSIM,
      nthin         = NTHIN,
      update_interval = 500
  )
  
  cate_draws_pathA <- fit_pathA$cate
  ate_draws_pathA  <- fit_pathA$ate
  
  cate_est_pathA   <- colMeans(cate_draws_pathA)
  cate_ci_pathA    <- apply(cate_draws_pathA, 2, quantile, probs = c(0.025, 0.975))
  
  rmse_pathA     <- sqrt(mean((cate_est_pathA - dgp$true_cate)^2))
  bias_pathA     <- mean(cate_est_pathA - dgp$true_cate)
  coverage_pathA <- mean(dgp$true_cate >= cate_ci_pathA[1, ] & dgp$true_cate <= cate_ci_pathA[2, ])
  cor_pathA      <- cor(cate_est_pathA, dgp$true_cate)

  # ----------------------------------------------------
  # Model 4: Tweedie BCF (Path B)
  # ----------------------------------------------------
  cat("\n  -> Fitting Tweedie BCF (Path B)...\n")
  fit_pathB <- countbcf_pathb(
      y             = dgp$y,
      z             = Z,
      x_control     = X,
      pihat         = pi_x,
      nburn         = NBURN,
      nsim          = NSIM,
      nthin         = NTHIN,
      update_interval = 500
  )
  
  # Potential outcomes: log-link mean parameters exp(mu_f)
  mu_f_B  <- fit_pathB$mu_f_post
  tau_f_B <- fit_pathB$tau_f_post
  
  cate_draws_pathB <- matrix(0, nrow = NSIM, ncol = N)
  for (s in 1:NSIM) {
    mu0_draw <- exp(2.0 * mu_f_B[s, ])
    mu1_draw <- exp(2.0 * (mu_f_B[s, ] + tau_f_B[s, ]))
    cate_draws_pathB[s, ] <- mu1_draw - mu0_draw
  }
  ate_draws_pathB <- rowMeans(cate_draws_pathB)
  
  cate_est_pathB  <- colMeans(cate_draws_pathB)
  cate_ci_pathB   <- apply(cate_draws_pathB, 2, quantile, probs = c(0.025, 0.975))
  
  rmse_pathB     <- sqrt(mean((cate_est_pathB - dgp$true_cate)^2))
  bias_pathB     <- mean(cate_est_pathB - dgp$true_cate)
  coverage_pathB <- mean(dgp$true_cate >= cate_ci_pathB[1, ] & dgp$true_cate <= cate_ci_pathB[2, ])
  cor_pathB      <- cor(cate_est_pathB, dgp$true_cate)
  
  # ----------------------------------------------------
  # Model 5: Joint Copula-BCF Path C (Selection model)
  # ----------------------------------------------------
  cat("\n  -> Fitting Joint Copula-BCF (Path C)...\n")
  fit_pathC <- pathc_bcf(
      y          = dgp$y,
      z          = Z,
      x_control  = X,
      pihat_sel  = pi_x,
      pihat_out  = NULL, # estimated automatically via SPA
      nburn      = NBURN,
      nsim       = NSIM,
      nthin      = NTHIN,
      update_interval = 500
  )
  
  # Re-transform Joint Copula-BCF draws to original scale
  # E[Y(z) | X] = exp(eta_c(z) + sigma^2 / 2) * Phi(eta_b(z) + beta)
  cate_draws_pathC <- matrix(0, nrow = NSIM, ncol = N)
  for (s in 1:NSIM) {
    eta_b0 <- fit_pathC$sel_con_post[s, ]
    eta_b1 <- fit_pathC$sel_con_post[s, ] + fit_pathC$sel_mod_post[s, ]
    
    eta_c0 <- fit_pathC$out_con_post[s, ]
    eta_c1 <- fit_pathC$out_con_post[s, ] + fit_pathC$out_mod_post[s, ]
    
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
  # Model 6: Gamma Hurdle BCF (Path D)
  # ----------------------------------------------------
  cat("\n  -> Fitting Gamma Hurdle BCF (Path D)...\n")
  fit_pathD_intensity <- pathd_gammabcf(
      y             = dgp$y,
      z             = Z,
      x_control     = X,
      pihat_pos     = NULL, # SPA automatically estimated
      nburn         = NBURN,
      nsim          = NSIM,
      thin          = NTHIN,
      update_interval = 500,
      return_trees  = TRUE
  )
  
  # Manually provide scale and shift parameters to fit_intensity control/moderate forests
  fit_pathD_intensity$control_fit$scale <- 1.0
  fit_pathD_intensity$control_fit$shift <- 0.0
  fit_pathD_intensity$moderate_fit$scale <- 1.0
  fit_pathD_intensity$moderate_fit$shift <- 0.0
  
  # Predict for all units using get_forest_fit, negating to obtain Conventional Log-Mean Scale
  X_c_all <- cbind(X, pi_x)
  mu_f_all  <- -get_forest_fit(fit_pathD_intensity$control_fit, X_c_all)   # nsim x n
  tau_f_all <- -get_forest_fit(fit_pathD_intensity$moderate_fit, X) # nsim x n
  
  # Re-use Probit Hurdle draws from Path A
  mu_b_all  <- fit_pathA$mu_b
  tau_b_all <- fit_pathA$tau_b
  
  p0_hurdle <- pnorm(mu_b_all)
  p1_hurdle <- pnorm(mu_b_all + tau_b_all)
  
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
  
  # Pack results
  data.frame(
    DGP = dgp$name,
    True_ATE = dgp$true_ate,
    Metric = c("Est ATE Mean", "Est ATE SD", "CATE RMSE", "CATE Abs Bias", "CATE 95% Coverage", "CATE Correlation"),
    BCF_Linear = c(mean(ate_draws_linear), sd(ate_draws_linear), rmse_linear, abs(bias_linear), coverage_linear, cor_linear),
    BCF_Log = c(mean(ate_draws_log), sd(ate_draws_log), rmse_log, abs(bias_log), coverage_log, cor_log),
    ZIC_BCF_PathA = c(mean(ate_draws_pathA), sd(ate_draws_pathA), rmse_pathA, abs(bias_pathA), coverage_pathA, cor_pathA),
    Tweedie_PathB = c(mean(ate_draws_pathB), sd(ate_draws_pathB), rmse_pathB, abs(bias_pathB), coverage_pathB, cor_pathB),
    Joint_Copula_PathC = c(mean(ate_draws_pathC), sd(ate_draws_pathC), rmse_pathC, abs(bias_pathC), coverage_pathC, cor_pathC),
    Gamma_Hurdle_PathD = c(mean(ate_draws_pathD), sd(ate_draws_pathD), rmse_pathD, abs(bias_pathD), coverage_pathD, cor_pathD)
  )
}

# Run all 3 DGPs
dgp_a <- generate_dgp_lognormal()
dgp_b <- generate_dgp_gaussian()
dgp_c <- generate_dgp_tweedie()

res_a <- run_eval(dgp_a)
res_b <- run_eval(dgp_b)
res_c <- run_eval(dgp_c)

# Combine and save results
all_results <- rbind(res_a, res_b, res_c)
write.csv(all_results, file.path(RESULTS_DIR, "full_simulation_results.csv"), row.names = FALSE)
cat(sprintf("\n[SUCCESS] Full comparative simulation results saved to: %s\n", file.path(RESULTS_DIR, "full_simulation_results.csv")))

# Print results nicely
print(all_results, digits = 4)
