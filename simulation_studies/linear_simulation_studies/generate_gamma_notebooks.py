"""
Generate self-contained Google Colab notebooks for the two parametric Gamma
benchmark models (Oganisian et al. 2019, appendix C): the Gamma hurdle model
and the naive Gamma +.01 model. Mirrors generate_benchmark_notebooks.py (the
EDP+ZI / DPglm generator); both Gamma models are shipped inside the zicbcf
package (functions gamma_hurdle() and gamma_plus01()), so each notebook only
needs devtools::install_github("hugogobato/zicbcf").

The DGP generators and calc_cate_metrics() are imported from
generate_benchmark_notebooks.py so they are byte-identical to every other
notebook and the benchmark numbers are directly comparable.

Unlike EDP, these GLM benchmarks are cheap (a couple of penalized-IRLS fits per
dataset, ~0.05-0.25 s each), so the full design is feasible to run on Colab.

Output layout (continuing the existing folder1..32 organisation):
  Gamma hurdle : folder33 standard | 34 n1000 | 35 n100 | 36 n250 | 37-40 zi lvl 1/2/4/5
  Gamma +.01   : folder41 standard | 42 n1000 | 43 n100 | 44 n250 | 45-48 zi lvl 1/2/4/5
"""
from generate_benchmark_notebooks import (
    create_notebook, install_cell, calc_metrics_r_code, dgps, zi_levels,
)

# gamma_hurdle has a hurdle component (P(Y>0) contrast); gamma_plus01 does not.
models = {
    "gamma_hurdle": {"name": "Gamma Hurdle (Oganisian et al. 2019)",
                     "prefix": "GammaHurdle", "code": "Gamma-Hurdle", "hurdle": True},
    "gamma_plus01": {"name": "Gamma +.01 (Oganisian et al. 2019)",
                     "prefix": "GammaP01", "code": "Gamma+.01", "hurdle": False},
}

CATE_FIELDS = ["CATE_RMSE", "CATE_Abs_Bias", "CATE_Coverage",
               "CATE_Correlation", "CATE_CI_Length", "Est_ATE"]
HURDLE_FIELDS = ["Hurdle_RMSE", "Hurdle_Abs_Bias", "Hurdle_Coverage",
                 "Hurdle_Correlation", "Hurdle_CI_Length", "Est_Hurdle_ATE"]
CATE_EXPR = ["m$rmse", "abs(m$bias)", "m$coverage", "m$correlation", "m$ci_length", "m$est_ate_mean"]
HURDLE_EXPR = ["mh$rmse", "abs(mh$bias)", "mh$coverage", "mh$correlation", "mh$ci_length", "mh$est_ate_mean"]


def _cols(pfx, fields, exprs):
    return [f"      {pfx}_{f} = {e}" for f, e in zip(fields, exprs)]


def fit_call(model_key):
    return (f"    {model_key}(y = d$y, z = d$z, x = d$x,\n"
            f"          {' ' * len(model_key)}nburn = NBURN, nsim = NSIM, nthin = NTHIN)")


def main_cell_src(model_key, dgp, c_shift_expr):
    m = models[model_key]
    pfx, name, has_hurdle = m["prefix"], m["name"], m["hurdle"]

    # OK branch: CATE always from m; hurdle from mh (if any) else NA.
    ok_cols = _cols(pfx, CATE_FIELDS, CATE_EXPR)
    if has_hurdle:
        ok_cols += _cols(pfx, HURDLE_FIELDS, HURDLE_EXPR)
        ok_hurdle_calc = (
            "    hurdle_cate_draws <- fit$p1 - fit$p0\n"
            "    mh <- calc_cate_metrics(hurdle_cate_draws, d$true_hurdle_cate, rowMeans(hurdle_cate_draws))\n"
        )
    else:
        ok_cols += _cols(pfx, HURDLE_FIELDS, ["NA"] * len(HURDLE_FIELDS))
        ok_hurdle_calc = ""
    ok_cols = ",\n".join(ok_cols)

    # FAILED branch: every metric column NA.
    na_cols = ",\n".join(_cols(pfx, CATE_FIELDS + HURDLE_FIELDS, ["NA"] * 12))

    return f"""cat("=== Starting {name} on {dgp['name']} ===\\n")
for (s in 1:N_SIM) {{
  cat(sprintf("[Seed %d/%d] Generating and fitting...\\n", s, N_SIM))
  d <- {dgp['func']}(N, P, seed = s, c_shift = {c_shift_expr})

  fit <- tryCatch(
{fit_call(model_key)},
    error = function(e) {{ message(sprintf("  seed %d FAILED: %s", s, conditionMessage(e))); NULL }}
  )

  if (is.null(fit)) {{
    # Record the failed seed transparently.
    df_res <- data.frame(
      DGP = "{dgp['name']}", Seed = s, Status = "FAILED",
      True_ATE = d$true_ate, True_Hurdle_ATE = d$true_hurdle_ate,
{na_cols},
      stringsAsFactors = FALSE)
  }} else {{
    m  <- calc_cate_metrics(fit$cate, d$true_cate, fit$ate)
{ok_hurdle_calc}    df_res <- data.frame(
      DGP = "{dgp['name']}", Seed = s, Status = "OK",
      True_ATE = d$true_ate, True_Hurdle_ATE = d$true_hurdle_ate,
{ok_cols},
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
    hurdle_note = ("CATE and Hurdle metrics" if model["hurdle"]
                   else "CATE metrics (Hurdle columns are NA: no hurdle component)")
    return [
        {"cell_type": "markdown", "id": "cell-title", "metadata": {}, "source": [
            f"# {title}\n",
            f"- **Model**: {model['name']} (`{model['code']}`) — benchmark, shipped in `zicbcf`\n",
            f"- **DGP**: {dgp['name']}\n",
            f"- **Sample Size (N)**: {n_val}\n",
            f"- **Zero-Inflation Intercept (c_shift)**: {c_shift_expr}\n",
            "- **MCMC**: NBURN = 1000, NSIM = 1000\n",
            "- **Simulations**: 100 seeds\n",
            f"- **Output**: CSV with {hurdle_note}; `Status` flags any failed seed.\n",
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
            dgps[dgp_key]["code"], "\n", calc_metrics_r_code,
        ]},
        {"cell_type": "code", "execution_count": None, "id": "cell-main", "metadata": {}, "source": [
            main_cell_src(model_key, dgp, c_shift_expr),
        ]},
    ]


model_start_folder = {"gamma_hurdle": 33, "gamma_plus01": 41}


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
    print(f"\nDone: {total} Gamma benchmark notebooks generated.")


if __name__ == "__main__":
    main()
