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

# Fit Gamma Intensity BCF directly using pathd_gammabcf
fit <- pathd_gammabcf(
  y = y, 
  z = z, 
  x_control = X, 
  nburn = 500, 
  nsim = 500,
  return_trees = TRUE
)

# Manually provide scale and shift parameters to fit control/moderate forests
fit$control_fit$scale <- 1.0
fit$control_fit$shift <- 0.0
fit$moderate_fit$scale <- 1.0
fit$moderate_fit$shift <- 0.0

active_idx <- fit$active_idx

# Test 1: predict using only X (this is what the old script did)
cat("--- Old Prediction Method (X only) ---\n")
mu_f_old <- -get_forest_fit(fit$control_fit, X[active_idx, , drop=FALSE])
cat("Old Predicted mu_f first 5 active obs means:\n")
print(colMeans(mu_f_old[, 1:5]))

# Test 2: predict using cbind(X, pihat_pos) for control, and X for moderate
cat("\n--- New Prediction Method (cbind(X, pihat_pos) for control) ---\n")
X_c_pos <- cbind(X[active_idx, , drop=FALSE], fit$pihat_pos)
mu_f_new <- -get_forest_fit(fit$control_fit, X_c_pos)
cat("New Predicted mu_f first 5 active obs means:\n")
print(colMeans(mu_f_new[, 1:5]))

cat("\nC++ actual fit$mu_f_post first 5 active obs means:\n")
print(colMeans(fit$mu_f_post[, 1:5]))
