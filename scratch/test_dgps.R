library(countbcf, lib.loc = "local_lib")

set.seed(42)
N     <- 1000
P     <- 5

X <- matrix(rnorm(N * P), N, P)
colnames(X) <- paste0("X", 1:P)
pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
Z    <- rbinom(N, 1, pi_x)

test_dgp_lognormal <- function(c_shift) {
  p_hurdle_0   <- pnorm(c_shift + 0.5 * X[, 1] - 0.3 * X[, 3])
  p_hurdle_1   <- pnorm(c_shift + 0.5 * X[, 1] - 0.3 * X[, 3] + 0.4 + 0.2 * X[, 1])
  p_hurdle_obs <- ifelse(Z == 1, p_hurdle_1, p_hurdle_0)
  I <- rbinom(N, 1, p_hurdle_obs)
  mean(I == 0)
}

test_dgp_tweedie <- function(c_shift) {
  log_mu0 <- 1.2 + c_shift + 0.8 * X[, 1] - 0.4 * X[, 3]
  log_mu1 <- 1.2 + c_shift + 0.8 * X[, 1] - 0.4 * X[, 3] + 0.6 + 0.3 * X[, 1]
  mu_true  <- ifelse(Z == 1, exp(log_mu1), exp(log_mu0))
  phi_true <- 1.5
  lambda_true <- 2 * sqrt(mu_true) / phi_true
  N_latent    <- rpois(N, lambda_true)
  mean(N_latent == 0)
}

cat("Log-Normal Zeros vs c_shift:\n")
for (c in c(-1.5, -1.0, -0.5, 0.0, 0.2, 0.5, 1.0, 1.5)) {
  cat(sprintf("c = %4.1f: %5.1f%%\n", c, 100 * test_dgp_lognormal(c)))
}

cat("\nTweedie Zeros vs c_shift:\n")
for (c in c(-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0)) {
  cat(sprintf("c = %4.1f: %5.1f%%\n", c, 100 * test_dgp_tweedie(c)))
}
