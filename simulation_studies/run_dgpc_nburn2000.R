################################################################################
##  DGP C (Tweedie Compound) Simulation Study with nburn=2000
##
##  Fits 5 models (excluding BCF-Log) for DGP C under nburn = 2000, nsim = 1000.
################################################################################

library(countbcf, lib.loc = "local_lib")

set.seed(42)

## ---- Global specs -----------------------------------------------------------
N     <- 1000
P     <- 5
NBURN <- 2000
NSIM  <- 1000
NTHIN <- 1

# Generate standard normal covariates
X <- matrix(rnorm(N * P), N, P)
colnames(X) <- paste0("X", 1:P)

# Propensity Score & Treatment Assignment (Confounded)
pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
Z    <- rbinom(N, 1, pi_x)

# DGP C: Tweedie Compound Poisson-Gamma DGP (p=1.5)
log_mu0 <- 1.2 + 0.8 * X[, 1] - 0.4 * X[, 3]
log_mu1 <- 1.2 + 0.8 * X[, 1] - 0.4 * X[, 3] + 0.6 + 0.3 * X[, 1]

mu0_true  <- exp(log_mu0)
mu1_true  <- exp(log_mu1)
true_cate <- mu1_true - mu0_true
true_ate  <- mean(true_cate)

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

dgp <- list(y = Y, true_cate = true_cate, true_ate = true_ate, name = "Tweedie Compound")

cat(sprintf("=== Evaluating DGP C: Tweedie Compound (True ATE = %.4f) with nburn=%d ===\n", true_ate, NBURN))

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
    update_interval = 9999
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
# Model 2: ZIC-BCF Path A (Two-part BCF with SPA)
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
    update_interval = 9999
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
# Model 3: Tweedie BCF (Path B)
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
    update_interval = 9999
)

# Potential outcomes: log-link mean parameters exp(2.0 * mu_f)
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
# Model 4: Joint Copula-BCF Path C (Selection model)
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
    update_interval = 9999
)

# Re-transform Joint Copula-BCF draws to original scale
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
# Model 5: Gamma Hurdle BCF (Path D)
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
    update_interval = 9999,
    return_trees  = TRUE
)

# Manually provide scale and shift parameters
fit_pathD_intensity$control_fit$scale <- 1.0
fit_pathD_intensity$control_fit$shift <- 0.0
fit_pathD_intensity$moderate_fit$scale <- 1.0
fit_pathD_intensity$moderate_fit$shift <- 0.0

mu_f_all  <- -get_forest_fit(fit_pathD_intensity$control_fit, cbind(X, pi_x))   # nsim x n
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

# Pack and print results nicely
out_df <- data.frame(
  Metric = c("Est ATE Mean", "Est ATE SD", "CATE RMSE", "CATE Abs Bias", "CATE 95% Coverage", "CATE Correlation"),
  BCF_Linear = c(mean(ate_draws_linear), sd(ate_draws_linear), rmse_linear, abs(bias_linear), coverage_linear, cor_linear),
  ZIC_BCF_PathA = c(mean(ate_draws_pathA), sd(ate_draws_pathA), rmse_pathA, abs(bias_pathA), coverage_pathA, cor_pathA),
  Tweedie_PathB = c(mean(ate_draws_pathB), sd(ate_draws_pathB), rmse_pathB, abs(bias_pathB), coverage_pathB, cor_pathB),
  Joint_Copula_PathC = c(mean(ate_draws_pathC), sd(ate_draws_pathC), rmse_pathC, abs(bias_pathC), coverage_pathC, cor_pathC),
  Gamma_Hurdle_PathD = c(mean(ate_draws_pathD), sd(ate_draws_pathD), rmse_pathD, abs(bias_pathD), coverage_pathD, cor_pathD)
)

print(out_df, digits = 4)
write.csv(out_df, "simulation_studies/results/dgpc_nburn2000_results.csv", row.names = FALSE)
cat("\n[SUCCESS] DGP C nburn=2000 results saved to simulation_studies/results/dgpc_nburn2000_results.csv\n")
