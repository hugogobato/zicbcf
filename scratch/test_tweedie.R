library(countbcf)

# Set seed
set.seed(42)

# Generate synthetic Tweedie data (p = 1.5)
n <- 500
p <- 5
X <- matrix(rnorm(n * p), nrow = n, ncol = p)

# Propensity score and treatment assignment
sigmoid <- function(z) 1 / (1 + exp(-z))
pi_x <- sigmoid(0.5 * X[, 1] - 0.4 * X[, 2])
Z <- rbinom(n, 1, pi_x)

# True mean functions (Tweedie mean parameter mu_i)
log_mu_0 <- 1.5 + 0.5 * X[, 1]
log_mu_1 <- log_mu_0 + 0.5 - 0.3 * X[, 2]

log_mu <- ifelse(Z == 1, log_mu_1, log_mu_0)
mu <- exp(log_mu)

# Tweedie p = 1.5 compound Poisson-exponential simulation
# True dispersion
phi_true <- 1.2

# Generate Tweedie outcomes
lambda <- 2 * sqrt(mu) / phi_true
rate <- 2 / (phi_true * sqrt(mu))

N_events <- rpois(n, lambda)
Y <- sapply(1:n, function(i) {
  if (N_events[i] == 0) {
    return(0)
  } else {
    return(sum(rexp(N_events[i], rate = rate[i])))
  }
})

cat("Realized fraction of zeros:", mean(Y == 0), "\n")
cat("Mean outcome (all):", mean(Y), "\n")
cat("Max outcome:", max(Y), "\n")

# Calculate true potential outcomes and true CATE
mu0_true <- exp(log_mu_0)
mu1_true <- exp(log_mu_1)
cate_true <- mu1_true - mu0_true
ate_true <- mean(cate_true)

cat("True ATE:", ate_true, "\n")

# Fit countbcf_pathb model with larger nburn and nsim to allow convergence
cat("Fitting countbcf_pathb...\n")
fit <- countbcf_pathb(
  y = Y,
  z = Z,
  x_control = X,
  pihat = pi_x,
  nburn = 2000,
  nsim = 1000,
  nthin = 1,
  update_interval = 250,
  ntree_control = 50,
  ntree_moderate = 20
)

# Extract predictions
mu_f_post <- fit$mu_f_post
tau_f_post <- fit$tau_f_post

# Correct scaling: trees predict 0.5 * log(mu), so we must multiply by 2 before exponentiating
mu0_post <- exp(2 * mu_f_post)
mu1_post <- exp(2 * (mu_f_post + tau_f_post))

cate_post <- mu1_post - mu0_post

mu0_hat <- colMeans(mu0_post)
mu1_hat <- colMeans(mu1_post)
cate_hat <- colMeans(cate_post)

cat("Estimated ATE (using mu_f_post & tau_f_post with factor 2):", mean(cate_hat), "\n")
cat("RMSE of CATE:", sqrt(mean((cate_hat - cate_true)^2)), "\n")
cat("Correlation between true and estimated CATE:", cor(cate_hat, cate_true), "\n")

# Let's also print Estimated ATE from fit$yhat directly
# fit$yhat has shape nsim x n, so observed_mu_hat should use colMeans
observed_mu_hat <- colMeans(fit$yhat)
observed_mu_true <- mu
cat("RMSE of observed mu:", sqrt(mean((observed_mu_hat - observed_mu_true)^2)), "\n")

# Let's print the posterior mean of phi
phi_post <- fit$phi
cat("True phi:", phi_true, "\n")
cat("Posterior mean of phi:", mean(phi_post), "\n")
cat("Posterior SD of phi:", sd(phi_post), "\n")
