################################################################################
##  Nonlinear versions of the three semicontinuous DGPs (A / B / C)
##
##  Motivation
##  ----------
##  In the main simulation study the Gamma-Hurdle benchmark is extremely
##  competitive. The conjecture is that this is an artefact of the DGPs being
##  *linear on the link scale*: every conditional-mean component in the original
##  generators is an affine function of X,
##      hurdle : pnorm(c_shift + 0.5*X1 - 0.3*X3   [+ treatment])
##      cont   : (log) 1.5 + 0.8*X2 + 0.4*X4       [+ treatment]
##  and the Gamma-Hurdle benchmark fits a fully treatment-interacted *linear*
##  design (1, z, X, z*X) with the matching links (logit + Gamma-log). With a
##  linear truth that design is (near-)correctly specified, so the parametric
##  GLM can essentially recover the true potential-outcome means -- and hence
##  the ATE and even the CATE -- leaving little room for a forest to win.
##
##  These nonlinear DGPs replace the affine predictors with functions the GLM
##  design CANNOT represent:
##    * covariate-by-covariate cross terms        (X1*X2, X2*X4, X1*X3)
##    * smooth univariate nonlinearities          (sin(c*X))
##    * centered quadratics                        (X^2 - 1)
##    * a treatment effect that is itself a cross term (z * X_i * X_j)
##  The Gamma-Hurdle design has only main effects + z + z*X, so it is genuinely
##  misspecified here, whereas BCF/forest models can adapt to the structure.
##
##  All nonlinear terms are mean-centered (sin of a symmetric variable, X_i*X_j
##  with i!=j, X^2-1) so the marginal zero proportion and outcome scale stay
##  close to the original "standard" configuration; the comparison therefore
##  isolates the effect of *functional form*, not of zero-inflation or scale.
##
##  The true-CATE / true-ATE formulas are unchanged in form -- they are computed
##  from the same (now nonlinear) predictors -- so the targets remain exact.
################################################################################

## ---- Shared nonlinear building blocks (mean ~ 0 for X ~ N(0,1)) --------------
## Coefficients are calibrated so the marginal zero proportion and outcome scale
## stay in the same ballpark as the original linear DGPs (A: ATE ~1.2, ~37% zeros;
## C: ~18% zeros). Treatment effects use BOUNDED nonlinearities (sin/tanh) so the
## CATE cannot blow up multiplicatively through the exp() link.
.nl_hurdle  <- function(X) 0.50 * sin(1.3 * X[, 1]) - 0.30 * X[, 3] +
                            0.35 * X[, 1] * X[, 2] - 0.20 * (X[, 3]^2 - 1)
.nl_tau_h   <- function(X) 0.20 + 0.20 * X[, 1] * X[, 2]                 # heterogeneous hurdle effect

.nl_cont    <- function(X) 1.50 + 0.80 * sin(1.2 * X[, 2]) + 0.40 * X[, 4] +
                            0.35 * X[, 2] * X[, 4] - 0.25 * (X[, 1]^2 - 1)
.nl_tau_c   <- function(X) 0.25 - 0.15 * sin(1.3 * X[, 2]) + 0.12 * tanh(X[, 1] * X[, 2])  # bounded intensity effect

## ============================================================================
##  DGP A (nonlinear): probit hurdle + Log-Normal continuous intensity
## ============================================================================
generate_dgp_a_nl <- function(n, p, seed, c_shift = 0.2) {
  set.seed(seed * 1000 + 42)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
  Z    <- rbinom(n, 1, pi_x)

  # Nonlinear hurdle (probit), halved-style heterogeneous treatment effect
  p_hurdle_0   <- pnorm(c_shift + .nl_hurdle(X))
  p_hurdle_1   <- pnorm(c_shift + .nl_hurdle(X) + .nl_tau_h(X))
  p_hurdle_obs <- ifelse(Z == 1, p_hurdle_1, p_hurdle_0)
  I <- rbinom(n, 1, p_hurdle_obs)

  # Nonlinear continuous intensity (Log-Normal)
  mu_c_0     <- .nl_cont(X)
  mu_c_1     <- .nl_cont(X) + .nl_tau_c(X)
  sigma_true <- 0.5

  y_pos_0   <- exp(mu_c_0 + rnorm(n, 0, sigma_true))
  y_pos_1   <- exp(mu_c_1 + rnorm(n, 0, sigma_true))
  y_pos_obs <- ifelse(Z == 1, y_pos_1, y_pos_0)
  Y <- I * y_pos_obs

  true_mu0  <- p_hurdle_0 * exp(mu_c_0 + 0.5 * sigma_true^2)
  true_mu1  <- p_hurdle_1 * exp(mu_c_1 + 0.5 * sigma_true^2)
  true_cate <- true_mu1 - true_mu0
  true_hurdle_cate <- p_hurdle_1 - p_hurdle_0

  list(y = Y, z = Z, x = X, pihat = pi_x,
       true_cate = true_cate, true_ate = mean(true_cate),
       true_hurdle_cate = true_hurdle_cate, true_hurdle_ate = mean(true_hurdle_cate))
}

## ============================================================================
##  DGP B (nonlinear): probit hurdle + Gamma continuous intensity (shape = 2)
## ============================================================================
generate_dgp_b_nl <- function(n, p, seed, c_shift = 0.2) {
  set.seed(seed * 1000 + 42)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
  Z    <- rbinom(n, 1, pi_x)

  p_hurdle_0   <- pnorm(c_shift + .nl_hurdle(X))
  p_hurdle_1   <- pnorm(c_shift + .nl_hurdle(X) + .nl_tau_h(X))
  p_hurdle_obs <- ifelse(Z == 1, p_hurdle_1, p_hurdle_0)
  I <- rbinom(n, 1, p_hurdle_obs)

  # Nonlinear log-mean for the Gamma intensity
  log_mu_c_0 <- .nl_cont(X)
  log_mu_c_1 <- .nl_cont(X) + .nl_tau_c(X)
  mu_c_0 <- exp(log_mu_c_0)
  mu_c_1 <- exp(log_mu_c_1)

  alpha   <- 2.0
  scale_0 <- mu_c_0 / alpha
  scale_1 <- mu_c_1 / alpha

  y_pos_0   <- rgamma(n, shape = alpha, scale = scale_0)
  y_pos_1   <- rgamma(n, shape = alpha, scale = scale_1)
  y_pos_obs <- ifelse(Z == 1, y_pos_1, y_pos_0)
  Y <- I * y_pos_obs

  true_mu0  <- p_hurdle_0 * mu_c_0
  true_mu1  <- p_hurdle_1 * mu_c_1
  true_cate <- true_mu1 - true_mu0
  true_hurdle_cate <- p_hurdle_1 - p_hurdle_0

  list(y = Y, z = Z, x = X, pihat = pi_x,
       true_cate = true_cate, true_ate = mean(true_cate),
       true_hurdle_cate = true_hurdle_cate, true_hurdle_ate = mean(true_hurdle_cate))
}

## ============================================================================
##  DGP C (nonlinear): compound Poisson-Gamma (Tweedie) with nonlinear log-mean
## ============================================================================
generate_dgp_c_nl <- function(n, p, seed, c_shift = 0.0) {
  set.seed(seed * 1000 + 42)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
  Z    <- rbinom(n, 1, pi_x)

  # Nonlinear log-mean (cross term X1*X3 + sin + centered quadratic)
  g_c   <- 1.2 + 0.80 * sin(1.2 * X[, 1]) - 0.40 * X[, 3] +
           0.35 * X[, 1] * X[, 3] - 0.25 * (X[, 2]^2 - 1)
  tau_c <- 0.30 - 0.15 * sin(1.3 * X[, 1]) + 0.12 * tanh(X[, 1] * X[, 3])  # bounded

  log_mu0 <- c_shift + g_c
  log_mu1 <- c_shift + g_c + tau_c

  mu0_true  <- exp(log_mu0)
  mu1_true  <- exp(log_mu1)
  true_cate <- mu1_true - mu0_true

  mu_true  <- ifelse(Z == 1, mu1_true, mu0_true)
  phi_true <- 1.5

  lambda0_true <- 2 * sqrt(mu0_true) / phi_true
  lambda1_true <- 2 * sqrt(mu1_true) / phi_true
  lambda_true  <- ifelse(Z == 1, lambda1_true, lambda0_true)

  N_latent   <- rpois(n, lambda_true)
  gamma_true <- 0.5 * phi_true * sqrt(mu_true)

  Y <- rep(0, n)
  for (i in 1:n) {
    if (N_latent[i] > 0) Y[i] <- rgamma(1, shape = N_latent[i], scale = gamma_true[i])
  }

  p0_hurdle_true   <- 1 - exp(-lambda0_true)
  p1_hurdle_true   <- 1 - exp(-lambda1_true)
  true_hurdle_cate <- p1_hurdle_true - p0_hurdle_true

  list(y = Y, z = Z, x = X, pihat = pi_x,
       true_cate = true_cate, true_ate = mean(true_cate),
       true_hurdle_cate = true_hurdle_cate, true_hurdle_ate = mean(true_hurdle_cate))
}

## ---- Shared CATE metric routine (identical to the main study) ----------------
calc_cate_metrics <- function(cate_draws, true_c, ate_draws) {
  cate_est <- colMeans(cate_draws)
  cate_ci  <- apply(cate_draws, 2, quantile, probs = c(0.025, 0.975))

  rmse <- sqrt(mean((cate_est - true_c)^2))
  bias <- mean(cate_est - true_c)
  coverage <- mean(true_c >= cate_ci[1, ] & true_c <= cate_ci[2, ])
  correlation <- cor(cate_est, true_c)
  if (is.na(correlation)) correlation <- 0.0
  ci_length <- mean(cate_ci[2, ] - cate_ci[1, ])
  est_ate_mean <- mean(ate_draws)

  list(rmse = rmse, bias = bias, coverage = coverage,
       correlation = correlation, ci_length = ci_length, est_ate_mean = est_ate_mean)
}

## Registry mirroring the main study's `dgps` (name / generator / standard c_shift).
## NOTE on DGP C: the nonlinear log-mean raises C's baseline, so the predictor's
## natural intercept (c_shift = 0) yields only ~10% zeros -- both too low for a
## "standard" cell and nearly tied with ZI level 4. We therefore set C's standard
## c_shift = -2.072, giving ~40% zeros (harmonized with A/B's ~40% standard). The
## ZI grid for C is correspondingly 85/60/40/11/3% (levels 1..5).
nl_dgps <- list(
  dgp_a = list(name = "DGP A (nonlinear): Log-Normal Hurdle",       func = generate_dgp_a_nl, c_shift = 0.2),
  dgp_b = list(name = "DGP B (nonlinear): Gamma Hurdle",            func = generate_dgp_b_nl, c_shift = 0.2),
  dgp_c = list(name = "DGP C (nonlinear): Tweedie Semicontinuous",  func = generate_dgp_c_nl, c_shift = -2.072)
)
