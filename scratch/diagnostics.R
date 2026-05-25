# Diagnostics for Path A with True Active Propensity Score
library(countbcf, lib.loc = "local_lib")

set.seed(42)

n <- 1000
p <- 5

X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)

# 1. Propensity Score & Treatment Assignment
pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
z <- rbinom(n, 1, pi_x)

# 2. Binary Participation Hurdle
p_hurdle_0 <- pnorm(0.2 + 0.5 * X[, 1] - 0.3 * X[, 3])
p_hurdle_1 <- pnorm(0.2 + 0.5 * X[, 1] - 0.3 * X[, 3] + 0.4 + 0.2 * X[, 1])
p_hurdle_obs <- ifelse(z == 1, p_hurdle_1, p_hurdle_0)

I <- rbinom(n, 1, p_hurdle_obs)
active_idx <- which(I == 1)

# 3. Continuous Intensity Component (Log-Normal)
mu_c_0 <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4]
mu_c_1 <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4] + 0.5 - 0.3 * X[, 2]
sigma_true <- 0.5

y_pos_0 <- exp(mu_c_0 + rnorm(n, 0, sigma_true))
y_pos_1 <- exp(mu_c_1 + rnorm(n, 0, sigma_true))
y_pos_obs <- ifelse(z == 1, y_pos_1, y_pos_0)

# 4. Final Semicontinuous Outcome
y <- I * y_pos_obs

# True active propensity score
p_I_z1 <- p_hurdle_1
p_I_z0 <- p_hurdle_0
pi_x_active_true_full <- (p_I_z1 * pi_x) / (p_I_z1 * pi_x + p_I_z0 * (1 - pi_x))
pi_x_active_true <- pi_x_active_true_full[active_idx]

# True response-scale values
true_mu0 <- p_hurdle_0 * exp(mu_c_0 + 0.5 * sigma_true^2)
true_mu1 <- p_hurdle_1 * exp(mu_c_1 + 0.5 * sigma_true^2)
true_cate <- true_mu1 - true_mu0
true_ate <- mean(true_cate)

cat("Fitting ZIC-BCF Path A with TRUE Active Propensity Score...\n")
fit <- zicbcf_pathA(
    y = y, z = z,
    x_control = X,
    pihat = pi_x,
    pihat_active = pi_x_active_true,
    nburn = 200, nsim = 500,
    update_interval = 100
)

cat("\n=== Component-level diagnostics ===\n")
# Hurdle Probit part
p0_est <- colMeans(pnorm(fit$mu_b))
p1_est <- colMeans(pnorm(fit$mu_b + fit$tau_b))
cat("Hurdle p0: True Mean =", mean(p_hurdle_0), "| Est Mean =", mean(p0_est), "\n")
cat("Hurdle p1: True Mean =", mean(p_hurdle_1), "| Est Mean =", mean(p1_est), "\n")

# Log-intensity part
mu_c_0_est <- colMeans(fit$mu_c)
mu_c_1_est <- colMeans(fit$mu_c + fit$tau_c)
cat("Log-intensity mu_c_0: True Mean =", mean(mu_c_0), "| Est Mean =", mean(mu_c_0_est), "\n")
cat("Log-intensity mu_c_1: True Mean =", mean(mu_c_1), "| Est Mean =", mean(mu_c_1_est), "\n")
cat("Sigma_c: True =", sigma_true, "| Est Mean =", mean(fit$sigma_c), "\n")

# Re-transformation details
y0_plus_est <- colMeans(exp(fit$mu_c + 0.5 * fit$sigma_c^2))
y1_plus_est <- colMeans(exp(fit$mu_c + fit$tau_c + 0.5 * fit$sigma_c^2))
cat("Response-scale y0_plus: True Mean =", mean(exp(mu_c_0 + 0.5*sigma_true^2)), "| Est Mean =", mean(y0_plus_est), "\n")
cat("Response-scale y1_plus: True Mean =", mean(exp(mu_c_1 + 0.5*sigma_true^2)), "| Est Mean =", mean(y1_plus_est), "\n")

# Full scale potential outcomes
mu0_est <- colMeans(fit$mu0)
mu1_est <- colMeans(fit$mu1)
cat("Response mu0: True Mean =", mean(true_mu0), "| Est Mean =", mean(mu0_est), "\n")
cat("Response mu1: True Mean =", mean(true_mu1), "| Est Mean =", mean(mu1_est), "\n")
cat("Response ATE: True =", true_ate, "| Est Mean =", mean(fit$ate), "\n")
