library(countbcf, lib.loc = "local_lib")

# Set seed
set.seed(42)

# Generate synthetic Tweedie data matching DGP C (Scenario 14)
N <- 1000
P <- 5
X <- matrix(rnorm(N * P), nrow = N, ncol = P)
colnames(X) <- paste0("X", 1:P)

# Propensity score and treatment assignment
pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
Z <- rbinom(N, 1, pi_x)

# True log-mean parameter
log_mu0 <- 1.2 + 0.8 * X[, 1] - 0.4 * X[, 3]
log_mu1 <- 1.2 + 0.8 * X[, 1] - 0.4 * X[, 3] + 0.6 + 0.3 * X[, 1]

mu0_true <- exp(log_mu0)
mu1_true <- exp(log_mu1)
true_cate <- mu1_true - mu0_true
true_ate <- mean(true_cate)

mu_true <- ifelse(Z == 1, mu1_true, mu0_true)
phi_true <- 1.5

# Simulation via latent Poisson and continuous Gamma sizes
lambda_true <- 2 * sqrt(mu_true) / phi_true
N_latent <- rpois(N, lambda_true)
gamma_true <- 0.5 * phi_true * sqrt(mu_true)

Y <- rep(0, N)
for (i in 1:N) {
  if (N_latent[i] > 0) {
    Y[i] <- rgamma(1, shape = N_latent[i], scale = gamma_true[i])
  }
}

cat("Realized fraction of zeros:", mean(Y == 0), "\n")
cat("Mean outcome (all):", mean(Y), "\n")
cat("True ATE:", true_ate, "\n")

# Fit Tweedie BCF (Path B) with long burn-in (8000 iterations)
cat("Fitting countbcf_pathb with nburn=8000 and nsim=2000...\n")
fit <- countbcf_pathb(
  y = Y,
  z = Z,
  x_control = X,
  pihat = pi_x,
  nburn = 8000,
  nsim = 2000,
  nthin = 1,
  update_interval = 1000
)

# Extract fits and correct potential outcomes scaling (factor of 2.0)
mu_f_post <- fit$mu_f_post
tau_f_post <- fit$tau_f_post

mu0_post <- exp(2.0 * mu_f_post)
mu1_post <- exp(2.0 * (mu_f_post + tau_f_post))

cate_post <- mu1_post - mu0_post
cate_hat <- colMeans(cate_post)
ate_hat <- rowMeans(cate_post)

cat("\n--- Tweedie BCF (Path B) Results ---\n")
cat("Estimated ATE (Mean):", mean(ate_hat), "\n")
cat("Estimated ATE (SD):", sd(ate_hat), "\n")
cat("Absolute Bias of ATE:", abs(mean(ate_hat) - true_ate), "\n")
cat("RMSE of CATE:", sqrt(mean((cate_hat - true_cate)^2)), "\n")
cat("Correlation between true and estimated CATE:", cor(cate_hat, true_cate), "\n")

# Coverage
cate_ci <- apply(cate_post, 2, quantile, probs = c(0.025, 0.975))
coverage <- mean(true_cate >= cate_ci[1, ] & true_cate <= cate_ci[2, ])
cat("CATE 95% Coverage:", coverage * 100, "%\n")

# Posterior of phi
cat("True phi:", phi_true, "\n")
cat("Posterior mean of phi:", mean(fit$phi), "\n")
cat("Posterior SD of phi:", sd(fit$phi), "\n")
