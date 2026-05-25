library(countbcf, lib.loc = "local_lib")
set.seed(42)
N <- 1000
P <- 5
X <- matrix(rnorm(N * P), N, P)
colnames(X) <- paste0("X", 1:P)
pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
Z    <- rbinom(N, 1, pi_x)

# Log-Normal Hurdle DGP (c=0.2)
p_hurdle_0   <- pnorm(0.2 + 0.5 * X[, 1] - 0.3 * X[, 3])
p_hurdle_1   <- pnorm(0.2 + 0.5 * X[, 1] - 0.3 * X[, 3] + 0.4 + 0.2 * X[, 1])
p_hurdle_obs <- ifelse(Z == 1, p_hurdle_1, p_hurdle_0)
I <- rbinom(N, 1, p_hurdle_obs)
mu_c_0     <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4]
mu_c_1     <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4] + 0.5 - 0.3 * X[, 2]
sigma_true <- 0.5
y_pos_0   <- exp(mu_c_0 + rnorm(N, 0, sigma_true))
y_pos_1   <- exp(mu_c_1 + rnorm(N, 0, sigma_true))
Y <- I * ifelse(Z == 1, y_pos_1, y_pos_0)

cat(sprintf("Active proportion: %.1f%%\n", 100 * mean(Y > 0)))

NBURN <- 100
NSIM  <- 200
NTHIN <- 1

methods <- list(
  "BCF-Linear" = function() {
    bcf_continuous_linear(
      y          = Y,
      z          = Z,
      x_control  = X,
      x_moderate = X,
      zhat       = pi_x,
      nburn      = NBURN,
      nsim       = NSIM,
      nthin      = NTHIN,
      update_interval = 500
    )
  },
  "BCF-Log" = function() {
    bcf_continuous_linear(
      y          = log(Y + 1),
      z          = Z,
      x_control  = X,
      x_moderate = X,
      zhat       = pi_x,
      nburn      = NBURN,
      nsim       = NSIM,
      nthin      = NTHIN,
      update_interval = 500
    )
  },
  "ZIC-BCF-PathA" = function() {
    zicbcf_pathA(
      y             = Y,
      z             = Z,
      x_control     = X,
      pihat         = pi_x,
      pihat_active  = NULL,
      nburn         = NBURN,
      nsim          = NSIM,
      nthin         = NTHIN,
      update_interval = 500
    )
  },
  "Tweedie-PathB" = function() {
    countbcf_pathb(
      y             = Y,
      z             = Z,
      x_control     = X,
      pihat         = pi_x,
      nburn         = NBURN,
      nsim          = NSIM,
      nthin         = NTHIN,
      update_interval = 500
    )
  },
  "Joint-Copula-PathC" = function() {
    pathc_bcf(
      y          = Y,
      z          = Z,
      x_control  = X,
      pihat_sel  = pi_x,
      pihat_out  = NULL,
      nburn      = NBURN,
      nsim       = NSIM,
      nthin      = NTHIN,
      update_interval = 500
    )
  },
  "Gamma-Hurdle-PathD" = function() {
    # Probit Hurdle is part of path A, but the Gamma intensity stage is path D
    pathd_gammabcf(
      y             = Y,
      z             = Z,
      x_control     = X,
      pihat_pos     = NULL,
      nburn         = NBURN,
      nsim          = NSIM,
      thin          = NTHIN,
      update_interval = 500,
      return_trees  = TRUE
    )
  }
)

for (m in names(methods)) {
  cat(sprintf("\nTiming method: %s...\n", m))
  t0 <- Sys.time()
  res <- methods[[m]]()
  print(Sys.time() - t0)
}
