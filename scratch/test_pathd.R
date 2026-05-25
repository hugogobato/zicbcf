library(countbcf)

set.seed(42)
n <- 1000

# 1. Generate baseline covariates
X <- matrix(rnorm(n * 5), n, 5)

# 2. Propensity score and treatment assignment
pi_score <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
z <- rbinom(n, 1, pi_score)

# 3. Binary participation hurdle (I = 1 if positive)
eta_b <- 0.2 + 0.5 * X[, 1] - 0.3 * X[, 3] + z * (0.4 + 0.2 * X[, 1])
prob_hurdle <- pnorm(eta_b)
I <- rbinom(n, 1, prob_hurdle)

# 4. Continuous intensity component (Gamma BCF)
kappa_true <- 4.0
log_lambda <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4] + z * (0.5 - 0.3 * X[, 2])
lambda <- exp(log_lambda)

# Positive continuous outcome conditional on participation
y_pos <- rgamma(n, shape = kappa_true, scale = lambda / kappa_true)

# Final semicontinuous outcome (Y = 0 if hurdle not passed)
y <- ifelse(I == 1, y_pos, 0.0)

# ==========================================
# FIT PATH D: GAMMA HURDLE BCF VIA GIG
# ==========================================

cat("Fitting intensity log-linear Gamma BCF...\n")
fit_intensity <- pathd_gammabcf(
  y = y, 
  z = z, 
  x_control = X, 
  pihat_pos = NULL, # SPA automatically estimated internally
  nburn = 500, 
  nsim = 500
)

# ==========================================
# CAUSAL INFERENCE & ESTIMAND EVALUATION
# ==========================================

# 1. Recover continuous component means
# Control potential outcome mean: lambda_0
# Treated potential outcome mean: lambda_1
lambda_0_draws <- fit_intensity$lambda_0_post
lambda_1_draws <- fit_intensity$lambda_1_post

cat("\nPath D MCMC completed successfully!\n")
cat("Estimated shape parameter kappa_c (mean posterior draw):", mean(fit_intensity$kappa), "\n")
cat("True shape parameter kappa_c was:", kappa_true, "\n")
cat("Summary of mu_f_post (first 5 obs means):", colMeans(fit_intensity$mu_f_post[, 1:5, drop=FALSE]), "\n")
cat("Summary of tau_f_post (first 5 obs means):", colMeans(fit_intensity$tau_f_post[, 1:5, drop=FALSE]), "\n")
