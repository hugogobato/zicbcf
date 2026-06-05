"""
Generate self-contained Google Colab notebooks for the NONLINEAR DGPs, mirroring
the structure of simulation_studies/colab_notebooks.

Models (5 notebook types, EDP+ZI excluded -- it is cubic-per-sweep / tens of
hours per notebook):
    bcf            -> BCF-Linear              (results_bcf_*)
    zicbcf_smear   -> ZIC-BCF-Smear           (results_zicbcf_smear_*)
    dpglm          -> DPglm                   (results_dpglm_*)
    gamma_hurdle   -> Gamma Hurdle            (results_gamma_hurdle_*)
    gamma_plus01   -> Gamma +.01              (results_gamma_plus01_*)

Configurations per DGP (8): standard (N=500) + N-sensitivity (1000/100/250) +
ZI-sensitivity (levels 1/2/4/5). MCMC budget = 1000 burn-in + 1000 saved draws,
100 seeds -- identical to the linear simulation studies.

Design notes
------------
* The DGP cell embeds nonlinear_dgps.R *verbatim* (single source of truth with
  run_nonlinear_comparison.R), so the notebooks fit exactly the same nonlinear
  DGPs we validated locally, and calc_cate_metrics() is byte-identical to the
  main study.
* The DPglm / Gamma main-code cells are produced by re-using the EXACT
  main_cell_src() builders from the existing simulation_studies generators, so
  the output columns (Status + *_CATE_* / *_Hurdle_*) match the linear study.
* ZI-level c_shifts were recalibrated (see calibrate_zi_cshifts below / README)
  so each nonlinear DGP hits the same target zero proportions as the linear
  study at ZI levels 1/2/4/5. The standard config (level 3) keeps the default
  c_shift (0.2 for A/B, 0.0 for C), matching run_nonlinear_comparison.R.

Output: nonlinear_simulation_studies/colab_notebooks/folderNN/<notebook>.ipynb
  bcf 1-8 | zicbcf_smear 9-16 | dpglm 17-24 | gamma_hurdle 25-32 | gamma_plus01 33-40
  (per model: 0 standard | 1 n1000 | 2 n100 | 3 n250 | 4-7 zi lvl 1/2/4/5)
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SIM  = os.path.join(REPO, "simulation_studies")
NB_BASE = os.path.join(HERE, "colab_notebooks")

# Reuse the existing builders (importing these modules has no side effects; their
# generation loops are guarded by `if __name__ == "__main__"`).
sys.path.insert(0, SIM)
import generate_benchmark_notebooks as gbn   # noqa: E402  (dpglm main_cell_src, install_cell)
import generate_gamma_notebooks as ggn       # noqa: E402  (gamma_hurdle / gamma_plus01 main_cell_src)

install_cell = gbn.install_cell

# Full nonlinear DGP source -> embedded as the cell-dgp of every notebook.
with open(os.path.join(HERE, "nonlinear_dgps.R"), encoding="utf-8") as f:
    NONLINEAR_DGP_R = f.read()

# --- Nonlinear DGP registry (names/func/standard-c_shift match nonlinear_dgps.R) -
# DGP C's standard c_shift is -2.072 (~40% zeros): the nonlinear log-mean raises
# C's baseline, so c_shift=0 gives only ~10% zeros. See nonlinear_dgps.R.
dgps = {
    "dgp_a": {"name": "DGP A (nonlinear): Log-Normal Hurdle",      "func": "generate_dgp_a_nl", "c_shift": "0.2"},
    "dgp_b": {"name": "DGP B (nonlinear): Gamma Hurdle",           "func": "generate_dgp_b_nl", "c_shift": "0.2"},
    "dgp_c": {"name": "DGP C (nonlinear): Tweedie Semicontinuous", "func": "generate_dgp_c_nl", "c_shift": "-2.072"},
}

# --- Recalibrated ZI-level c_shifts (realized zero prop ~ target per level) -----
# A/B follow the linear study's targets (85/62/[40 std]/18/7%). DGP C's grid is
# shifted up to span a proper range with a ~40% standard: 85/60/[40 std]/11/3%
#   level 1 (~85/85/85%), 2 (~62/61/60%), 4 (~18/18/11%), 5 (~7/5/3%)
zi_levels = {
    1: {"dgp_a": "-1.367", "dgp_b": "-1.352", "dgp_c": "-5.648"},
    2: {"dgp_a": "-0.465", "dgp_b": "-0.426", "dgp_c": "-3.315"},
    4: {"dgp_a": "1.005",  "dgp_b": "1.023",  "dgp_c": "-0.163"},
    5: {"dgp_a": "1.783",  "dgp_b": "1.978",  "dgp_c": "1.131"},
}

# --- Model registry: display name, code tag, start folder ----------------------
models = {
    "bcf":          {"name": "BCF-Linear",        "code": "BCF-Linear",    "start": 1},
    "zicbcf_smear": {"name": "ZIC-BCF-Smear",     "code": "ZIC-BCF-Smear", "start": 9},
    "dpglm":        {"name": "DPglm (Oganisian et al. 2021)", "code": "DPglm", "start": 17},
    "gamma_hurdle": {"name": "Gamma Hurdle (Oganisian et al. 2019)", "code": "Gamma-Hurdle", "start": 25},
    "gamma_plus01": {"name": "Gamma +.01 (Oganisian et al. 2019)",   "code": "Gamma+.01",    "start": 33},
}


def create_notebook(folder, filename, cells):
    d = os.path.join(NB_BASE, folder)
    os.makedirs(d, exist_ok=True)
    nb = {
        "nbformat": 4, "nbformat_minor": 5,
        "metadata": {"kernelspec": {"display_name": "R", "language": "R", "name": "ir"},
                     "language_info": {"name": "R"}},
        "cells": cells,
    }
    with open(os.path.join(d, filename), "w", encoding="utf-8") as f:
        json.dump(nb, f, indent=1)


# --- Main-code cell for the two BCF-family notebooks --------------------------
def bcf_linear_main(dgp, c_shift_expr):
    return f"""cat("=== Starting BCF-Linear Simulation ===\\n")
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
cat("\\n=== Finished BCF-Linear Run ===\\n")"""


def smear_main(dgp, c_shift_expr):
    return f"""cat("=== Starting ZIC-BCF-Smear Simulation ===\\n")
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


def main_code(model_key, dgp, c_shift_expr):
    """Dispatch to the right main-code-cell builder for each model."""
    if model_key == "bcf":
        return bcf_linear_main(dgp, c_shift_expr)
    if model_key == "zicbcf_smear":
        return smear_main(dgp, c_shift_expr)
    if model_key == "dpglm":
        return gbn.main_cell_src("dpglm", dgp, c_shift_expr)
    if model_key in ("gamma_hurdle", "gamma_plus01"):
        return ggn.main_cell_src(model_key, dgp, c_shift_expr)
    raise ValueError(model_key)


def build_cells(model_key, dgp_key, n_val, c_shift_expr, lvl=None):
    dgp = dgps[dgp_key]
    m = models[model_key]
    title = f"Colab Simulation (Nonlinear DGP): {m['name']} on {dgp['name']}"
    # ZI-sensitivity cells get a level suffix so they never clobber the standard
    # (level 3) N=500 file; standard + N-sensitivity keep the plain name.
    out_csv = (f"results_{model_key}_{dgp_key}_N{n_val}.csv" if lvl is None
               else f"results_{model_key}_{dgp_key}_N{n_val}_lvl_{lvl}.csv")
    return [
        {"cell_type": "markdown", "id": "cell-title", "metadata": {}, "source": [
            f"# {title}\n",
            f"- **Model**: {m['name']} (`{m['code']}`)\n",
            f"- **DGP**: {dgp['name']} (nonlinear conditional means)\n",
            f"- **Sample Size (N)**: {n_val}\n",
            f"- **Zero-Inflation Intercept (c_shift)**: {c_shift_expr}\n",
            "- **MCMC**: NBURN = 1000, NSIM = 1000\n",
            "- **Simulations**: 100 seeds\n",
            "- **Output**: CSV with CATE and Hurdle metrics; `Status` flags any failed seed (benchmarks).\n",
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
            f"OUT_CSV <- \"{out_csv}\"\n",
            "if (file.exists(OUT_CSV)) file.remove(OUT_CSV)",
        ]},
        {"cell_type": "code", "execution_count": None, "id": "cell-dgp", "metadata": {}, "source": [
            NONLINEAR_DGP_R,
        ]},
        {"cell_type": "code", "execution_count": None, "id": "cell-main", "metadata": {}, "source": [
            main_code(model_key, dgp, c_shift_expr),
        ]},
    ]


# Per-model config order: (folder_offset, kind, params)
CONFIGS = [
    (0, "standard", {}),
    (1, "n", {"n_val": 1000}),
    (2, "n", {"n_val": 100}),
    (3, "n", {"n_val": 250}),
    (4, "zi", {"lvl": 1}),
    (5, "zi", {"lvl": 2}),
    (6, "zi", {"lvl": 4}),
    (7, "zi", {"lvl": 5}),
]


def file_stub(model_key, kind, n_val=None, lvl=None):
    if kind == "standard":
        return f"run_{model_key}_"
    if kind == "n":
        return f"run_{model_key}_sensitivity_n{n_val}_"
    return f"run_{model_key}_sensitivity_zi_lvl{lvl}_"


def main():
    total = 0
    for model_key, m in models.items():
        f0 = m["start"]
        for off, kind, prm in CONFIGS:
            folder = f"folder{f0 + off}"
            for dgp_key, dgp in dgps.items():
                lvl = None
                if kind == "standard":
                    n_val, c_shift = 500, dgp["c_shift"]
                    stub = file_stub(model_key, "standard")
                elif kind == "n":
                    n_val, c_shift = prm["n_val"], dgp["c_shift"]
                    stub = file_stub(model_key, "n", n_val=n_val)
                else:
                    lvl = prm["lvl"]
                    n_val, c_shift = 500, zi_levels[lvl][dgp_key]
                    stub = file_stub(model_key, "zi", lvl=lvl)
                fname = f"{stub}{dgp_key}.ipynb"
                create_notebook(folder, fname, build_cells(model_key, dgp_key, n_val, c_shift, lvl=lvl))
                total += 1
    print(f"Done: {total} nonlinear notebooks generated under {NB_BASE}")


if __name__ == "__main__":
    main()
