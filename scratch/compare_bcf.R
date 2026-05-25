# Compare ZIC-BCF continuous stage with standard bcf_continuous_linear
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

# Observed semicontinuous y
y <- I * y_pos_obs

# True active propensity score
p_I_z1 <- p_hurdle_1
p_I_z0 <- p_hurdle_0
pi_x_active_true_full <- (p_I_z1 * pi_x) / (p_I_z1 * pi_x + p_I_z0 * (1 - pi_x))
pi_x_active_true <- pi_x_active_true_full[active_idx]

# Let's fit standard continuous BCF on log(y_pos_obs[active_idx])
y_cont <- log(y[active_idx])
X_active <- X[active_idx, ]
z_active <- z[active_idx]

cat("Fitting standard bcf_continuous_linear on active subset...\n")
fit_bcf <- bcf_continuous_linear(
    y = y_cont,
    z = z_active,
    x_control = X_active,
    x_moderate = X_active,
    zhat = pi_x_active_true,
    nburn = 200, nsim = 500,
    update_interval = 100
)

# Recover mu0 and mu1 on log scale from standard BCF:
muy_active <- mean(y_cont)
mu_post <- muy_active + get_forest_fit(fit_bcf$control_fit, X_active)
tau_post <- get_forest_fit(fit_bcf$moderate_fit, X_active)

mu0_bcf <- colMeans(mu_post)
cate_bcf <- colMeans(tau_post)
mu1_bcf <- mu0_bcf + cate_bcf

cat("\n=== Standard BCF Results on Active Subset ===\n")
cat("True mu_c_0 mean:       ", mean(mu_c_0[active_idx]), "\n")
cat("Estimated mu_c_0 mean:  ", mean(mu0_bcf), "\n")
cat("True mu_c_1 mean:       ", mean(mu_c_1[active_idx]), "\n")
cat("Estimated mu_c_1 mean:  ", mean(mu1_bcf), "\n")
cat("True CATE (log) mean:   ", mean((mu_c_1 - mu_c_0)[active_idx]), "\n")
cat("Estimated CATE mean:    ", mean(cate_bcf), "\n")
cat("=========================================\n")
