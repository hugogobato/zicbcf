import json
import os

# Create folder for Colab notebooks
notebooks_dir = "/home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks"
os.makedirs(notebooks_dir, exist_ok=True)

def create_notebook(filename, cells):
    nb = {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {
                "display_name": "R",
                "language": "R",
                "name": "ir"
            },
            "language_info": {
                "name": "R"
            }
        },
        "cells": cells
    }
    filepath = os.path.join(notebooks_dir, filename)
    with open(filepath, 'w', encoding='utf-8') as f:
        json.dump(nb, f, indent=1)
    print(f"Created notebook: {filename}")

# Shared R helper libraries installation cell
install_cell = {
    "cell_type": "code",
    "execution_count": None,
    "id": "cell-install",
    "metadata": {},
    "outputs": [],
    "source": [
        "# Install devtools and the zicbcf package from GitHub\n",
        "install.packages(\"remotes\", repos=\"https://cloud.r-project.org/\")\n",
        "if (!require(\"devtools\")) {\n",
        "  install.packages(\"devtools\", repos=\"https://cloud.r-project.org/\")\n",
        "}\n",
        "devtools::install_github(\"hugogobato/zicbcf\")\n",
        "library(zicbcf)"
    ]
}

# R code for the DGP generators taking N and c_shift as parameters
dgp_a_code = """generate_dgp_a <- function(n, p, seed, c_shift = 0.2) {
  set.seed(seed * 1000 + 42)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  
  pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
  Z    <- rbinom(n, 1, pi_x)
  
  # Hurdle (halved treatment effect)
  p_hurdle_0   <- pnorm(c_shift + 0.5 * X[, 1] - 0.3 * X[, 3])
  p_hurdle_1   <- pnorm(c_shift + 0.5 * X[, 1] - 0.3 * X[, 3] + 0.2 + 0.1 * X[, 1])
  p_hurdle_obs <- ifelse(Z == 1, p_hurdle_1, p_hurdle_0)
  I <- rbinom(n, 1, p_hurdle_obs)
  
  # Continuous intensity Log-Normal (halved treatment effect)
  mu_c_0     <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4]
  mu_c_1     <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4] + 0.25 - 0.15 * X[, 2]
  sigma_true <- 0.5
  
  y_pos_0   <- exp(mu_c_0 + rnorm(n, 0, sigma_true))
  y_pos_1   <- exp(mu_c_1 + rnorm(n, 0, sigma_true))
  y_pos_obs <- ifelse(Z == 1, y_pos_1, y_pos_0)
  Y <- I * y_pos_obs
  
  # True potential outcomes & CATE
  true_mu0  <- p_hurdle_0 * exp(mu_c_0 + 0.5 * sigma_true^2)
  true_mu1  <- p_hurdle_1 * exp(mu_c_1 + 0.5 * sigma_true^2)
  true_cate <- true_mu1 - true_mu0
  
  true_hurdle_cate <- p_hurdle_1 - p_hurdle_0
  
  list(y = Y, z = Z, x = X, pihat = pi_x, true_cate = true_cate, true_ate = mean(true_cate), 
       true_hurdle_cate = true_hurdle_cate, true_hurdle_ate = mean(true_hurdle_cate))
}"""

dgp_b_code = """generate_dgp_b <- function(n, p, seed, c_shift = 0.2) {
  set.seed(seed * 1000 + 42)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  
  pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
  Z    <- rbinom(n, 1, pi_x)
  
  # Hurdle (similar treatment effect to A, halved)
  p_hurdle_0   <- pnorm(c_shift + 0.5 * X[, 1] - 0.3 * X[, 3])
  p_hurdle_1   <- pnorm(c_shift + 0.5 * X[, 1] - 0.3 * X[, 3] + 0.2 + 0.1 * X[, 1])
  p_hurdle_obs <- ifelse(Z == 1, p_hurdle_1, p_hurdle_0)
  I <- rbinom(n, 1, p_hurdle_obs)
  
  # Continuous intensity Gamma (right-skewed, zero-inflated, different from A & C)
  log_mu_c_0 <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4]
  log_mu_c_1 <- 1.5 + 0.8 * X[, 2] + 0.4 * X[, 4] + 0.25 - 0.15 * X[, 2]
  
  mu_c_0 <- exp(log_mu_c_0)
  mu_c_1 <- exp(log_mu_c_1)
  
  alpha <- 2.0 # shape parameter for right-skewness
  scale_0 <- mu_c_0 / alpha
  scale_1 <- mu_c_1 / alpha
  
  y_pos_0 <- rgamma(n, shape = alpha, scale = scale_0)
  y_pos_1 <- rgamma(n, shape = alpha, scale = scale_1)
  y_pos_obs <- ifelse(Z == 1, y_pos_1, y_pos_0)
  Y <- I * y_pos_obs
  
  # True potential outcomes & CATE
  true_mu0  <- p_hurdle_0 * mu_c_0
  true_mu1  <- p_hurdle_1 * mu_c_1
  true_cate <- true_mu1 - true_mu0
  
  true_hurdle_cate <- p_hurdle_1 - p_hurdle_0
  
  list(y = Y, z = Z, x = X, pihat = pi_x, true_cate = true_cate, true_ate = mean(true_cate), 
       true_hurdle_cate = true_hurdle_cate, true_hurdle_ate = mean(true_hurdle_cate))
}"""

dgp_c_code = """generate_dgp_c <- function(n, p, seed, c_shift = 0.0) {
  set.seed(seed * 1000 + 42)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  
  pi_x <- pnorm(-0.5 + 0.4 * X[, 1] + 0.3 * X[, 2]^2)
  Z    <- rbinom(n, 1, pi_x)
  
  # Log-mean parameter (halved treatment effect)
  log_mu0 <- 1.2 + c_shift + 0.8 * X[, 1] - 0.4 * X[, 3]
  log_mu1 <- 1.2 + c_shift + 0.8 * X[, 1] - 0.4 * X[, 3] + 0.3 + 0.15 * X[, 1]
  
  mu0_true  <- exp(log_mu0)
  mu1_true  <- exp(log_mu1)
  true_cate <- mu1_true - mu0_true
  
  mu_true  <- ifelse(Z == 1, mu1_true, mu0_true)
  phi_true <- 1.5
  
  lambda0_true <- 2 * sqrt(mu0_true) / phi_true
  lambda1_true <- 2 * sqrt(mu1_true) / phi_true
  
  lambda_true <- ifelse(Z == 1, lambda1_true, lambda0_true)
  N_latent    <- rpois(n, lambda_true)
  gamma_true  <- 0.5 * phi_true * sqrt(mu_true)
  
  Y <- rep(0, n)
  for (i in 1:n) {
    if (N_latent[i] > 0) {
      Y[i] <- rgamma(1, shape = N_latent[i], scale = gamma_true[i])
    }
  }
  
  p0_hurdle_true <- 1 - exp(-lambda0_true)
  p1_hurdle_true <- 1 - exp(-lambda1_true)
  true_hurdle_cate <- p1_hurdle_true - p0_hurdle_true
  
  list(y = Y, z = Z, x = X, pihat = pi_x, true_cate = true_cate, true_ate = mean(true_cate), 
       true_hurdle_cate = true_hurdle_cate, true_hurdle_ate = mean(true_hurdle_cate))
}"""

# Shared evaluation metric collector (precisely matching run_5seed_single.R)
calc_metrics_r_code = """calc_cate_metrics <- function(cate_draws, true_c, ate_draws) {
  cate_est <- colMeans(cate_draws)
  cate_ci  <- apply(cate_draws, 2, quantile, probs = c(0.025, 0.975))
  
  rmse <- sqrt(mean((cate_est - true_c)^2))
  bias <- mean(cate_est - true_c)
  coverage <- mean(true_c >= cate_ci[1, ] & true_c <= cate_ci[2, ])
  correlation <- cor(cate_est, true_c)
  if (is.na(correlation)) correlation <- 0.0
  ci_length <- mean(cate_ci[2, ] - cate_ci[1, ])
  est_ate_mean <- mean(ate_draws)
  
  list(rmse=rmse, bias=bias, coverage=coverage, correlation=correlation, ci_length=ci_length, est_ate_mean=est_ate_mean)
}"""

# ------------------------------------------------------------------------------
# 1. BUILD PRIMARY NOTEBOOK SPECIFICATIONS
# ------------------------------------------------------------------------------

dgps = {
    "dgp_a": {"name": "DGP A: Log-Normal Hurdle", "func": "generate_dgp_a", "c_shift": "0.2", "code": dgp_a_code},
    "dgp_b": {"name": "DGP B: Gamma Hurdle", "func": "generate_dgp_b", "c_shift": "0.2", "code": dgp_b_code},
    "dgp_c": {"name": "DGP C: Tweedie Semicontinuous", "func": "generate_dgp_c", "c_shift": "0.0", "code": dgp_c_code}
}

models = {
    "bcf": {"name": "Standard Continuous Gaussian BCF", "m_code": "BCF-Linear"},
    "zicbcf_smear": {"name": "ZIC-BCF-Smear", "m_code": "ZIC-BCF-Smear"}
}

def generate_notebook_cells(model_key, dgp_key, n_val, c_shift_expr):
    dgp = dgps[dgp_key]
    model = models[model_key]
    
    title = f"Colab Simulation: {model['name']} on {dgp['name']}"
    subtitle = f"- **Sample Size (N)**: {n_val}\n- **Zero-Inflation Intercept (c_shift)**: {c_shift_expr}\n"
    
    # R main execution logic
    main_execution = ""
    if model_key == "bcf":
        main_execution = f"""cat("=== Starting Standard BCF Simulation ===\\n")
for (s in 1:N_SIM) {{
  cat(sprintf("[Seed %d/%d] Generating and fitting...\\n", s, N_SIM))
  d <- {dgp['func']}(N, P, seed = s, c_shift = {c_shift_expr})
  
  fit <- bcf_continuous_linear(
    y          = d$y,
    z          = d$z,
    x_control  = d$x,
    x_moderate = d$x,
    zhat       = d$pihat,
    nburn      = NBURN,
    nsim       = NSIM,
    nthin      = NTHIN,
    update_interval = 99999
  )
  
  cate_draws <- get_forest_fit(fit$moderate_fit, d$x)
  ate_draws  <- rowMeans(cate_draws)
  m <- calc_cate_metrics(cate_draws, d$true_cate, ate_draws)
  
  df_res <- data.frame(
    DGP = "{dgp['name']}",
    Seed = s,
    True_ATE = d$true_ate,
    True_Hurdle_ATE = d$true_hurdle_ate,
    
    Linear_CATE_RMSE        = m$rmse,
    Linear_CATE_Abs_Bias    = abs(m$bias),
    Linear_CATE_Coverage    = m$coverage,
    Linear_CATE_Correlation = m$correlation,
    Linear_CATE_CI_Length   = m$ci_length,
    Linear_Est_ATE          = m$est_ate_mean,
    
    # Smear columns populated as NA for BCF-Linear notebooks
    Smear_CATE_RMSE         = NA,
    Smear_CATE_Abs_Bias     = NA,
    Smear_CATE_Coverage     = NA,
    Smear_CATE_Correlation  = NA,
    Smear_CATE_CI_Length    = NA,
    Smear_Est_ATE           = NA,
    
    Linear_Hurdle_RMSE        = NA,
    Linear_Hurdle_Abs_Bias    = NA,
    Linear_Hurdle_Coverage    = NA,
    Linear_Hurdle_Correlation = NA,
    Linear_Hurdle_CI_Length   = NA,
    Linear_Est_Hurdle_ATE     = NA,
    
    Smear_Hurdle_RMSE        = NA,
    Smear_Hurdle_Abs_Bias    = NA,
    Smear_Hurdle_Coverage    = NA,
    Smear_Hurdle_Correlation = NA,
    Smear_Hurdle_CI_Length   = NA,
    Smear_Est_Hurdle_ATE     = NA,
    stringsAsFactors = FALSE
  )
  
  write.table(df_res, OUT_CSV, sep=",", row.names=FALSE, col.names=!file.exists(OUT_CSV), append=TRUE)
  gc(verbose=FALSE)
}}
cat("\\n=== Finished Standard BCF Run ===\\n")"""
    else:
        main_execution = f"""cat("=== Starting ZIC-BCF-Smear Simulation ===\\n")
for (s in 1:N_SIM) {{
  cat(sprintf("[Seed %d/%d] Generating and fitting...\\n", s, N_SIM))
  d <- {dgp['func']}(N, P, seed = s, c_shift = {c_shift_expr})
  
  fit <- zicbcf_smear(
    y             = d$y,
    z             = d$z,
    x_control     = d$x,
    pihat         = d$pihat,
    nburn         = NBURN,
    nsim          = NSIM,
    update_interval = 99999
  )
  
  m <- calc_cate_metrics(fit$cate, d$true_cate, fit$ate)
  
  # Probit Hurdle stage draws
  p0_draws <- pnorm(fit$mu_b)
  p1_draws <- pnorm(fit$mu_b + fit$tau_b)
  hurdle_cate_draws <- p1_draws - p0_draws
  hurdle_ate_draws  <- rowMeans(hurdle_cate_draws)
  m_hurdle <- calc_cate_metrics(hurdle_cate_draws, d$true_hurdle_cate, hurdle_ate_draws)
  
  df_res <- data.frame(
    DGP = "{dgp['name']}",
    Seed = s,
    True_ATE = d$true_ate,
    True_Hurdle_ATE = d$true_hurdle_ate,
    
    Linear_CATE_RMSE         = NA,
    Linear_CATE_Abs_Bias     = NA,
    Linear_CATE_Coverage     = NA,
    Linear_CATE_Correlation  = NA,
    Linear_CATE_CI_Length    = NA,
    Linear_Est_ATE           = NA,
    
    Smear_CATE_RMSE        = m$rmse,
    Smear_CATE_Abs_Bias    = abs(m$bias),
    Smear_CATE_Coverage    = m$coverage,
    Smear_CATE_Correlation = m$correlation,
    Smear_CATE_CI_Length   = m$ci_length,
    Smear_Est_ATE          = m$est_ate_mean,
    
    Linear_Hurdle_RMSE        = NA,
    Linear_Hurdle_Abs_Bias    = NA,
    Linear_Hurdle_Coverage    = NA,
    Linear_Hurdle_Correlation = NA,
    Linear_Hurdle_CI_Length   = NA,
    Linear_Est_Hurdle_ATE     = NA,
    
    Smear_Hurdle_RMSE        = m_hurdle$rmse,
    Smear_Hurdle_Abs_Bias    = abs(m_hurdle$bias),
    Smear_Hurdle_Coverage    = m_hurdle$coverage,
    Smear_Hurdle_Correlation = m_hurdle$correlation,
    Smear_Hurdle_CI_Length   = m_hurdle$ci_length,
    Smear_Est_Hurdle_ATE     = m_hurdle$est_ate_mean,
    stringsAsFactors = FALSE
  )
  
  write.table(df_res, OUT_CSV, sep=",", row.names=FALSE, col.names=!file.exists(OUT_CSV), append=TRUE)
  gc(verbose=FALSE)
}}
cat("\\n=== Finished ZIC-BCF-Smear Run ===\\n")"""

    cells = [
        {
            "cell_type": "markdown",
            "id": "cell-title",
            "metadata": {},
            "source": [
                f"# {title}\n",
                f"- **Model**: {model['name']} (`{model['m_code']}`)\n",
                f"- **DGP**: {dgp['name']}\n",
                f"{subtitle}",
                "- **MCMC**: NBURN = 1000, NSIM = 1000\n",
                "- **Simulations**: 100 seeds\n",
                "- **Output**: CSV containing CATE and Hurdle metrics.\n"
            ]
        },
        install_cell,
        {
            "cell_type": "code",
            "execution_count": None,
            "id": "cell-params",
            "metadata": {},
            "source": [
                "N_SIM   <- 100L\n",
                f"N       <- {n_val}L\n",
                "P       <- 5L\n",
                "NBURN   <- 1000L\n",
                "NSIM    <- 1000L\n",
                "NTHIN   <- 1L\n",
                "\n",
                f"OUT_CSV <- \"results_{model_key}_{dgp_key}_N{n_val}.csv\"\n",
                "if (file.exists(OUT_CSV)) file.remove(OUT_CSV)"
            ]
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "id": "cell-dgp",
            "metadata": {},
            "source": [
                dgp['code'],
                "\n",
                calc_metrics_r_code
            ]
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "id": "cell-main",
            "metadata": {},
            "source": [
                main_execution
            ]
        }
    ]
    return cells

# ------------------------------------------------------------------------------
# 2. WRITE OUT ALL THE NOTEBOOKS
# ------------------------------------------------------------------------------

# 2.1 Standard runs: N = 500, c_shift = standard (0.2 for A/B, 0.0 for C)
for m_key in models:
    for dgp_key, dgp in dgps.items():
        filename = f"run_{m_key}_{dgp_key}.ipynb"
        cells = generate_notebook_cells(m_key, dgp_key, 500, dgp['c_shift'])
        create_notebook(filename, cells)

# 2.2 Sample Size (N) Sensitivity: N = (100, 250, 1000) (standard c_shifts)
for n_val in [100, 250, 1000]:
    for m_key in models:
        for dgp_key, dgp in dgps.items():
            filename = f"run_{m_key}_sensitivity_n{n_val}_{dgp_key}.ipynb"
            cells = generate_notebook_cells(m_key, dgp_key, n_val, dgp['c_shift'])
            create_notebook(filename, cells)

# 2.3 Zero-Inflation (ZI) Sensitivity: N = 500, c_shift = Level 1, 2, 4, 5
zi_levels = {
    1: {"dgp_a": "-1.5", "dgp_b": "-1.5", "dgp_c": "-3.5"},
    2: {"dgp_a": "-0.5", "dgp_b": "-0.5", "dgp_c": "-2.0"},
    4: {"dgp_a": "1.0",  "dgp_b": "1.0",  "dgp_c": "0.0"},
    5: {"dgp_a": "1.8",  "dgp_b": "1.8",  "dgp_c": "1.0"}
}

for lvl, shifts in zi_levels.items():
    for m_key in models:
        for dgp_key, dgp in dgps.items():
            filename = f"run_{m_key}_sensitivity_zi_lvl{lvl}_{dgp_key}.ipynb"
            cells = generate_notebook_cells(m_key, dgp_key, 500, shifts[dgp_key])
            create_notebook(filename, cells)

print("Programmatic generation complete.")
