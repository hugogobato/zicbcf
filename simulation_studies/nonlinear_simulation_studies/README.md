# Nonlinear-DGP Stress Test: ZIC-BCF-Smear vs. Gamma Hurdle

*Does the Gamma-Hurdle benchmark only look strong because the simulation DGPs are linear?*

Last updated: 2026-06-03.

---

## 1. The hypothesis

In the main simulation study the parametric **Gamma Hurdle** benchmark is
surprisingly competitive with the proposed **ZIC-BCF-Smear** — sometimes winning
on ATE, and even on CATE in the Tweedie DGP. The conjecture tested here is that
this is an **artefact of the DGPs being linear on the link scale**.

Every conditional-mean component in the original generators is an *affine*
function of the covariates:

```
hurdle (probit) :  c_shift + 0.5*X1 - 0.3*X3        [+ treatment]
continuous (log):  1.5 + 0.8*X2 + 0.4*X4            [+ treatment]
```

The Gamma-Hurdle benchmark fits a fully treatment-interacted **linear** design
`(1, z, X, z·X)` with the matching links (logistic hurdle + Gamma log-link). When
the truth is link-linear, that design is **(near-)correctly specified**, so the
GLM can essentially recover the true potential-outcome means — and hence the ATE
and even the CATE — leaving little room for a forest to win.

A natural prediction follows (and matches the intuition that *"the ATE is
recoverable but the CATE is not"*): the **ATE** is a global average and is
relatively forgiving of functional-form error, whereas the **CATE** (per-unit
heterogeneity) is not. So even a mild nonlinearity should hurt Gamma Hurdle's
CATE far more than its ATE.

## 2. What the nonlinear DGPs change

`nonlinear_dgps.R` keeps the three DGP *families* (A = log-normal hurdle,
B = Gamma hurdle, C = Tweedie) and their positive-part distributions, but
replaces the affine predictors with functions the GLM design **cannot
represent**:

- covariate-by-covariate **cross terms**: `X1*X2`, `X2*X4`, `X1*X3`
- smooth **univariate nonlinearities**: `sin(c·X)`
- **centered quadratics**: `X^2 - 1`
- a **treatment effect that is itself nonlinear/interacted** (`z` × `sin`/`tanh`
  of cross terms)

The Gamma-Hurdle design has only main effects + `z` + `z·X`, so none of these
terms are in its column space — it is genuinely misspecified. BCF/forests, by
contrast, can split on interactions and nonlinearities.

**Calibration.** All nonlinear terms are mean-centered (`sin` of a symmetric
variable, `Xi·Xj` with `i≠j`, `X²−1`), and the treatment effects use *bounded*
nonlinearities (`sin`, `tanh`) so the CATE cannot blow up multiplicatively
through the `exp()` link. The result keeps the outcome scale in the same ballpark
as the original linear DGPs. At the standard configuration all three DGPs sit at
≈ 40% zeros (A/B via the default `c_shift`; C via `c_shift = -2.072`, since its
nonlinear log-mean would otherwise give only ~10% zeros — see §7).
The comparison therefore isolates the effect of **functional form**, not of
zero-inflation or scale. The true-CATE/ATE formulas are unchanged in form (they
are computed from the same, now nonlinear, predictors), so the targets stay
exact.

## 3. Files

| File | Purpose |
|------|---------|
| `nonlinear_dgps.R` | Nonlinear `generate_dgp_{a,b,c}_nl()` generators + `calc_cate_metrics()` (identical metric routine to the main study) + the `nl_dgps` registry. |
| `run_nonlinear_comparison.R` | Quick local comparison (ZIC-BCF-Smear vs Gamma Hurdle); writes `results/`. |
| `results/nonlinear_results_long.csv` | One row per (DGP, seed, model) with all CATE + hurdle metrics. |
| `results/nonlinear_summary.csv` | Aggregated means/RMSE over seeds. |
| `generate_nonlinear_notebooks.py` | Generates the full Colab notebook grid under `colab_notebooks/` (see §7). |
| `run_nonlinear_notebooks_local.py` | Runs those notebooks locally (drops the GitHub-install cell, isolates output per folder); `--smoke` for a fast end-to-end check. |
| `colab_notebooks/folder1..40/` | The generated notebooks (5 models × 3 DGPs × 8 configs = 120). |

### Reproduce

```bash
cd nonlinear_simulation_studies
Rscript run_nonlinear_comparison.R 3     # 3 seeds/DGP (default); pass e.g. 20 for more
```

Standard configuration only: **N = 500**, default `c_shift`, MCMC **1000 burn-in
+ 1000 saved draws** — the identical budget to the linear simulation studies.
Each `zicbcf_smear` fit is ~30 s; `gamma_hurdle` is ~0.1 s, so 3 seeds × 3 DGPs
runs in ~5 minutes.

## 4. Results (3 seeds per DGP, standard configuration)

Mean over the 3 seeds (ATE RMSE = root-mean-square of the per-seed ATE error):

| DGP | Model | CATE RMSE | CATE Corr | CATE Cov | CATE CI Len | ATE RMSE |
| :-- | :-- | :--: | :--: | :--: | :--: | :--: |
| **A: Log-Normal** | Gamma Hurdle | 4.924 | **−0.388** | 66.6% | 6.750 | **0.393** |
| | **ZIC-BCF-Smear** | **1.662** | **+0.688** | 86.5% | 3.924 | 0.618 |
| **B: Gamma** | Gamma Hurdle | 4.318 | **−0.192** | 69.0% | 6.659 | **0.517** |
| | **ZIC-BCF-Smear** | **1.817** | **+0.651** | 85.3% | 4.551 | 0.892 |
| **C: Tweedie** | Gamma Hurdle | 0.311 | **−0.105** | 89.9% | 0.934 | 0.133 |
| | **ZIC-BCF-Smear** | **0.172** | **+0.163** | 94.3% | 0.639 | **0.084** |

*(DGP C's standard is now ~40% zeros — see §7; its absolute CATE/ATE scale is
smaller than A/B because the larger zero mass shrinks the response-scale effects.
Smear wins C on **every** metric, including ATE.)*

## 5. Interpretation — the hypothesis is confirmed

The ranking from the main (linear) study **flips** once the DGPs are nonlinear:

1. **CATE collapses for Gamma Hurdle.** Its CATE RMSE is ~1.8–3× worse than
   ZIC-BCF-Smear's in every DGP. The decisive signal is the **CATE correlation,
   which goes negative** in all three (A −0.39, B −0.19, C −0.11): the linear
   design does not merely lose precision, it ranks units' treatment effects in
   roughly the *wrong order*. ZIC-BCF-Smear keeps a positive correlation
   (0.69 / 0.65 / 0.16).

2. **Coverage degrades for Gamma Hurdle** (67–90%) because its posterior
   contracts around a misspecified mean — confidently wrong. ZIC-BCF-Smear stays
   closer to nominal (85–94%).

3. **The ATE is forgiving — exactly as predicted.** On the *global average*
   effect, Gamma Hurdle remains competitive (even better than Smear in A and B:
   0.39 vs 0.62, 0.52 vs 0.89), because functional-form errors largely average
   out over the population. Only in the Tweedie DGP (C) does Smear also win the
   ATE. This is the crux: **a link-linear GLM can still estimate the ATE under
   nonlinearity, but not the CATE.**

**Conclusion.** The Gamma-Hurdle benchmark's headline strength in the main study
is, to a large extent, a **linearity artefact**. It is a strong comparator *only*
when the conditional-mean surface is (close to) linear on the link scale and/or
when the ATE — not the CATE — is the estimand. As soon as the heterogeneity is a
genuinely nonlinear function of the covariates, the parametric design recovers
the wrong heterogeneity, and the nonparametric **ZIC-BCF-Smear** dominates on
CATE while remaining competitive on ATE.

## 6. Caveats

- **3 seeds is indicative, not definitive.** It is a stress test / proof of
  concept; the effect sizes here are large and consistent across all three DGPs
  and across CATE RMSE / correlation / coverage, but a full study (≥100 seeds,
  and the BCF-Linear / DPglm / Gamma +.01 comparators, N- and ZI-sensitivity
  grids) would be needed for paper-grade numbers. Re-run with a larger seed count
  via `Rscript run_nonlinear_comparison.R 100`.
- ZIC-BCF-Smear's coverage here (~85–92%) is mildly below nominal — the nonlinear
  DGPs are harder and 3 seeds are noisy — but it is far better calibrated than
  Gamma Hurdle's.
- The nonlinear coefficient choices are one reasonable instantiation; the
  qualitative conclusion (CATE collapse vs. ATE robustness for the linear GLM) is
  the robust takeaway, not the exact magnitudes.

## 7. Full Colab notebook grid (`colab_notebooks/`)

For a paper-grade study (100 seeds, all comparators, full sensitivity grid),
`generate_nonlinear_notebooks.py` produces the same notebook layout as
`simulation_studies/colab_notebooks`, but on the **nonlinear** DGPs:

- **5 models** (EDP+ZI excluded — cubic per sweep, tens of hours/notebook):
  `bcf` (BCF-Linear), `zicbcf_smear`, `dpglm`, `gamma_hurdle`, `gamma_plus01`.
- **3 DGPs** × **8 configs**: standard (N=500) + N-sensitivity (N=100/250/1000)
  + ZI-sensitivity (levels 1/2/4/5) → **120 notebooks**.
- MCMC = 1000 burn-in + 1000 saved draws, 100 seeds (identical to the linear
  study). Each is self-contained: `devtools::install_github("hugogobato/zicbcf")`
  then run top-to-bottom; the DGP cell embeds `nonlinear_dgps.R` verbatim.

Folder map (per model: `0` standard · `1/2/3` N=1000/100/250 · `4/5/6/7` ZI lvl 1/2/4/5):

| Model | Folders |
|-------|---------|
| BCF-Linear | `folder1`–`8` |
| ZIC-BCF-Smear | `folder9`–`16` |
| DPglm | `folder17`–`24` |
| Gamma Hurdle | `folder25`–`32` |
| Gamma +.01 | `folder33`–`40` |

Output CSVs follow `results_<model>_<dgp>_N<n>.csv` for the standard and
N-sensitivity cells, and `results_<model>_<dgp>_N500_lvl_<k>.csv` for the
ZI-sensitivity levels (k = 1/2/4/5). The level suffix makes every filename in the
grid unique, so all results can share one flat directory (the standard N=500 file
is reserved for level 3). The column **schema** matches the linear study, so
`simulation_studies/generate_report_elements.R` can aggregate them with only the
`RESULTS_DIR` / DGP-name / ZI-path strings adjusted.

### ZI-level zero proportions (recalibrated)

The ZI `c_shift`s were recalibrated (`uniroot` on a 200k-row draw) so each
nonlinear DGP hits the linear study's target zero proportions; the standard
config keeps the default `c_shift` (matching §4):

| Level | `c_shift` A / B / C | realized zeros A / B / C |
|------|------|------|
| 1 | −1.367 / −1.352 / −5.648 | 85% / 85% / 85% |
| 2 | −0.465 / −0.426 / −3.315 | 62% / 61% / 60% |
| **3 (standard)** | **0.2 / 0.2 / −2.072** | **40% / 40% / 40%** |
| 4 | 1.005 / 1.023 / −0.163 | 18% / 18% / 11% |
| 5 | 1.783 / 1.978 / 1.131 | 7% / 5% / 3% |

*(DGP C's grid was shifted up so its standard sits at ~40% zeros — matching A/B —
rather than the ~10% the default `c_shift = 0` would give: its old level 1/2 became
level 2/3 and a new ~85% level 1 was added. The grid is monotone in all three DGPs.)*

### Generate / run

```bash
cd nonlinear_simulation_studies
python3 generate_nonlinear_notebooks.py                      # (re)build the 120 notebooks
python3 run_nonlinear_notebooks_local.py --folders 9 --smoke # fast end-to-end check (N_SIM=1)
python3 run_nonlinear_notebooks_local.py --jobs 6            # full local run -> results_notebooks/
```

Smoke-tested end-to-end (1 seed, reduced MCMC) for one folder per model: all five
model types run cleanly and write the expected schema. The `gamma_hurdle`
notebook already reproduces the negative CATE correlation seen in §4.
