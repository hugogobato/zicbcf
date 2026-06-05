"""
Generate self-contained Google Colab notebooks for the EDP+ZI and DPglm
benchmark models (Kim et al. 2024), mirroring the existing BCF / ZIC-BCF-Smear
notebooks in colab_notebooks/. Both benchmarks are now shipped inside the
zicbcf package (functions edp_zi() and dpglm()), so each notebook only needs
devtools::install_github("hugogobato/zicbcf").

The DGP generators and calc_cate_metrics() are byte-identical to the ones used
by the existing notebooks so that the benchmark results are directly comparable.

Output layout (continuing the existing folder1..16 organisation):
  EDP+ZI : folder17 standard | 18 n1000 | 19 n100 | 20 n250 | 21-24 zi lvl 1/2/4/5
  DPglm  : folder25 standard | 26 n1000 | 27 n100 | 28 n250 | 29-32 zi lvl 1/2/4/5
"""
import json
import os

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "colab_notebooks")


def create_notebook(folder, filename, cells):
    d = os.path.join(BASE, folder)
    os.makedirs(d, exist_ok=True)
    nb = {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {"display_name": "R", "language": "R", "name": "ir"},
            "language_info": {"name": "R"},
        },
        "cells": cells,
    }
    with open(os.path.join(d, filename), "w", encoding="utf-8") as f:
        json.dump(nb, f, indent=1)
    print(f"Created {folder}/{filename}")


install_cell = {
    "cell_type": "code", "execution_count": None, "id": "cell-install", "metadata": {}, "outputs": [],
    "source": [
        "# Install the zicbcf package from GitHub (ships BCF, ZIC-BCF, and the\n",
        "# EDP+ZI / DPglm benchmark models with all their dependencies).\n",
        "install.packages(\"remotes\", repos=\"https://cloud.r-project.org/\")\n",
        "if (!require(\"devtools\")) {\n",
        "  install.packages(\"devtools\", repos=\"https://cloud.r-project.org/\")\n",
        "}\n",
        "devtools::install_github(\"hugogobato/zicbcf\")\n",
        "library(zicbcf)",
    ],
}

# --- DGP generators (identical to generate_colab_notebooks_final.py) ----------
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

dgps = {
    "dgp_a": {"name": "DGP A: Log-Normal Hurdle", "func": "generate_dgp_a", "c_shift": "0.2", "code": dgp_a_code},
    "dgp_b": {"name": "DGP B: Gamma Hurdle",       "func": "generate_dgp_b", "c_shift": "0.2", "code": dgp_b_code},
    "dgp_c": {"name": "DGP C: Tweedie Semicontinuous", "func": "generate_dgp_c", "c_shift": "0.0", "code": dgp_c_code},
}

models = {
    "edp":   {"name": "EDP+ZI (Kim et al. 2024)", "prefix": "EDP",   "code": "EDP-ZI"},
    "dpglm": {"name": "DPglm (Oganisian et al. 2021)", "prefix": "DPglm", "code": "DPglm"},
}


def fit_call(model_key):
    if model_key == "edp":
        return ("    edp_zi(y = d$y, z = d$z, x = d$x,\n"
                "           nburn = NBURN, nsim = NSIM, nthin = NTHIN, update_interval = 99999)")
    return ("    dpglm(y = d$y, z = d$z, x = d$x,\n"
            "          nburn = NBURN, nsim = NSIM, trim = NTHIN)")


def main_cell_src(model_key, dgp, c_shift_expr):
    pfx = models[model_key]["prefix"]
    name = models[model_key]["name"]
    # NA-filled metric column block, reused for failed seeds and present for all rows
    cols = (
        f"      {pfx}_CATE_RMSE        = m$rmse,\n"
        f"      {pfx}_CATE_Abs_Bias    = abs(m$bias),\n"
        f"      {pfx}_CATE_Coverage    = m$coverage,\n"
        f"      {pfx}_CATE_Correlation = m$correlation,\n"
        f"      {pfx}_CATE_CI_Length   = m$ci_length,\n"
        f"      {pfx}_Est_ATE          = m$est_ate_mean,\n"
        f"      {pfx}_Hurdle_RMSE        = mh$rmse,\n"
        f"      {pfx}_Hurdle_Abs_Bias    = abs(mh$bias),\n"
        f"      {pfx}_Hurdle_Coverage    = mh$coverage,\n"
        f"      {pfx}_Hurdle_Correlation = mh$correlation,\n"
        f"      {pfx}_Hurdle_CI_Length   = mh$ci_length,\n"
        f"      {pfx}_Est_Hurdle_ATE     = mh$est_ate_mean"
    )
    na_cols = "\n".join(
        f"      {pfx}_{c} = NA," if not c.endswith("Est_Hurdle_ATE") else f"      {pfx}_{c} = NA"
        for c in ["CATE_RMSE", "CATE_Abs_Bias", "CATE_Coverage", "CATE_Correlation", "CATE_CI_Length", "Est_ATE",
                  "Hurdle_RMSE", "Hurdle_Abs_Bias", "Hurdle_Coverage", "Hurdle_Correlation", "Hurdle_CI_Length", "Est_Hurdle_ATE"]
    )
    return f"""cat("=== Starting {name} on {dgp['name']} ===\\n")
for (s in 1:N_SIM) {{
  cat(sprintf("[Seed %d/%d] Generating and fitting...\\n", s, N_SIM))
  d <- {dgp['func']}(N, P, seed = s, c_shift = {c_shift_expr})

  fit <- tryCatch(
{fit_call(model_key)},
    error = function(e) {{ message(sprintf("  seed %d FAILED: %s", s, conditionMessage(e))); NULL }}
  )

  if (is.null(fit)) {{
    # Record the failed seed transparently (mirrors the authors' tryCatch-skip).
    df_res <- data.frame(
      DGP = "{dgp['name']}", Seed = s, Status = "FAILED",
      True_ATE = d$true_ate, True_Hurdle_ATE = d$true_hurdle_ate,
{na_cols},
      stringsAsFactors = FALSE)
  }} else {{
    m  <- calc_cate_metrics(fit$cate, d$true_cate, fit$ate)
    hurdle_cate_draws <- fit$p1 - fit$p0
    mh <- calc_cate_metrics(hurdle_cate_draws, d$true_hurdle_cate, rowMeans(hurdle_cate_draws))
    df_res <- data.frame(
      DGP = "{dgp['name']}", Seed = s, Status = "OK",
      True_ATE = d$true_ate, True_Hurdle_ATE = d$true_hurdle_ate,
{cols},
      stringsAsFactors = FALSE)
  }}

  write.table(df_res, OUT_CSV, sep=",", row.names=FALSE,
              col.names=!file.exists(OUT_CSV), append=TRUE)
  gc(verbose=FALSE)
}}
cat("\\n=== Finished {name} Run ===\\n")"""


def build_cells(model_key, dgp_key, n_val, c_shift_expr):
    dgp = dgps[dgp_key]
    model = models[model_key]
    title = f"Colab Simulation: {model['name']} on {dgp['name']}"
    return [
        {"cell_type": "markdown", "id": "cell-title", "metadata": {}, "source": [
            f"# {title}\n",
            f"- **Model**: {model['name']} (`{model['code']}`) — benchmark, shipped in `zicbcf`\n",
            f"- **DGP**: {dgp['name']}\n",
            f"- **Sample Size (N)**: {n_val}\n",
            f"- **Zero-Inflation Intercept (c_shift)**: {c_shift_expr}\n",
            "- **MCMC**: NBURN = 1000, NSIM = 1000\n",
            "- **Simulations**: 100 seeds\n",
            "- **Output**: CSV with CATE and Hurdle metrics; `Status` flags any failed seed.\n",
        ]},
        install_cell,
        {"cell_type": "code", "execution_count": None, "id": "cell-params", "metadata": {}, "source": [
            "N_SIM   <- 100L\n",
            f"N       <- {n_val}L\n",
            "P       <- 5L\n",
            "NBURN   <- 1000L\n",
            "NSIM    <- 1000L\n",
            "NTHIN   <- 1L\n",
            "\n",
            f"OUT_CSV <- \"results_{model_key}_{dgp_key}_N{n_val}.csv\"\n",
            "if (file.exists(OUT_CSV)) file.remove(OUT_CSV)",
        ]},
        {"cell_type": "code", "execution_count": None, "id": "cell-dgp", "metadata": {}, "source": [
            dgp["code"], "\n", calc_metrics_r_code,
        ]},
        {"cell_type": "code", "execution_count": None, "id": "cell-main", "metadata": {}, "source": [
            main_cell_src(model_key, dgp, c_shift_expr),
        ]},
    ]


# Folder plan: per model, 8 configs in the existing order
#   [standard, n1000, n100, n250, zi1, zi2, zi4, zi5]
zi_levels = {
    1: {"dgp_a": "-1.5", "dgp_b": "-1.5", "dgp_c": "-3.5"},
    2: {"dgp_a": "-0.5", "dgp_b": "-0.5", "dgp_c": "-2.0"},
    4: {"dgp_a": "1.0",  "dgp_b": "1.0",  "dgp_c": "0.0"},
    5: {"dgp_a": "1.8",  "dgp_b": "1.8",  "dgp_c": "1.0"},
}

model_start_folder = {"edp": 17, "dpglm": 25}


def file_stub(model_key, kind, n_val=None, lvl=None):
    if kind == "standard":
        return f"run_{model_key}_"
    if kind == "n":
        return f"run_{model_key}_sensitivity_n{n_val}_"
    return f"run_{model_key}_sensitivity_zi_lvl{lvl}_"


def main():
    total = 0
    for model_key in models:
        f0 = model_start_folder[model_key]
        # config list: (folder_offset, kind, params)
        configs = [
            (0, "standard", {}),
            (1, "n", {"n_val": 1000}),
            (2, "n", {"n_val": 100}),
            (3, "n", {"n_val": 250}),
            (4, "zi", {"lvl": 1}),
            (5, "zi", {"lvl": 2}),
            (6, "zi", {"lvl": 4}),
            (7, "zi", {"lvl": 5}),
        ]
        for off, kind, prm in configs:
            folder = f"folder{f0 + off}"
            for dgp_key, dgp in dgps.items():
                if kind == "standard":
                    n_val, c_shift = 500, dgp["c_shift"]
                    stub = file_stub(model_key, "standard")
                elif kind == "n":
                    n_val, c_shift = prm["n_val"], dgp["c_shift"]
                    stub = file_stub(model_key, "n", n_val=n_val)
                else:
                    n_val, c_shift = 500, zi_levels[prm["lvl"]][dgp_key]
                    stub = file_stub(model_key, "zi", lvl=prm["lvl"])
                fname = f"{stub}{dgp_key}.ipynb"
                create_notebook(folder, fname, build_cells(model_key, dgp_key, n_val, c_shift))
                total += 1
    print(f"\nDone: {total} benchmark notebooks generated under {BASE}")


if __name__ == "__main__":
    main()
