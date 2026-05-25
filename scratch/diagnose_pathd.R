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

# Active subset
active_idx <- which(y > 0)
y_pos_active <- y[active_idx]
z_pos_active <- z[active_idx]
X_pos_active <- X[active_idx, , drop=FALSE]

# Fit Gamma Intensity BCF directly using pathd_gammabcf
fit <- pathd_gammabcf(
  y = y, 
  z = z, 
  x_control = X, 
  nburn = 500, 
  nsim = 500
)

cat("--- Raw C++ and transformed output comparison ---\n")
cat("True log_lambda for first 5 active obs:\n")
print(log_lambda[active_idx[1:5]])

cat("\nfit$yhat_log (transformed log_lambda) first 5 active obs means:\n")
print(colMeans(fit$yhat_log[, 1:5]))

cat("\nfit$mu_f_post (transformed control log_lambda) first 5 active obs means:\n")
print(colMeans(fit$mu_f_post[, 1:5]))

cat("\nfit$tau_f_post (transformed treatment effect) first 5 active obs means:\n")
print(colMeans(fit$tau_f_post[, 1:5]))
