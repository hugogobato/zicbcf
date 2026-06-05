# Integration of the Gamma Hurdle and Gamma +.01 Benchmark Models

*Documentation of how the two parametric Gamma benchmark models from Oganisian
et al. (2019) were implemented and integrated into the `zicbcf` package and the
simulation study. Companion to `benchmark_models_EDP_DPglm.md`; written to
support the methods / reproducibility sections of the paper.*

Last updated: 2026-06-02.

---

## 1. Objective

To complement the Bayesian-nonparametric benchmarks (EDP+ZI, DPglm) with the two
**parametric** comparators used by Oganisian, Mitra & Roy (2019), *"A Bayesian
Nonparametric Model for Zero-Inflated Outcomes: Prediction, Clustering, and
Causal Estimation"* (arXiv:1810.09494v3), appendix C:

1. **Gamma hurdle** — a two-part model: a logistic regression for the
   probability of a positive outcome \(P(Y>0)\), and a Gamma regression with a
   log link for the positive part \(Y\mid Y>0\).
2. **Gamma +.01** — the naive trick of adding a small constant (\(0.01\)) to the
   zeros and fitting a single Gamma log-link regression to \(Y+0.01\).

Both are run on the **same** DGPs (A/B/C) and the same sample-size and
zero-inflation sensitivity grid as the proposed estimators, with directly
comparable outputs, and runnable on Colab through a single
`devtools::install_github("hugogobato/zicbcf")`.

## 2. What the paper specifies (and what it leaves open)

The paper states only that these models are "coded ... in Stan" with Bayesian
posterior draws (5000 after 5000 burn-in), used as **ATE** comparators
(Table 1 reports bias / coverage / interval width of the average effect). It
does not give the design or priors for the simulation fits, and it never reports
a per-subject CATE or a zero-probability contrast for these comparators.

Two modelling choices were therefore made explicit for our CATE-focused study,
and are documented here:

- **Design.** Both parts use a fully treatment-by-covariate interacted design
  \((1, z, X, z\!\cdot\!X)\). The main-effects-only specification a researcher
  might fit for an ATE cannot express a covariate-varying treatment effect; the
  interacted design lets the parametric model produce a genuine heterogeneous
  CATE (so the correlation/coverage metrics are a fair, non-strawman comparison).
  Response-scale effects then follow by **G-computation (standardization)**:
  \(E[Y\mid do(a),x]=\mathrm{expit}(\eta_a)\cdot\exp(\xi_a)\) for the hurdle, and
  \(\exp(\xi_a)\) for +.01.
- **Fit / priors.** Each part is fit as a **weakly-informative Bayesian GLM** —
  see §3.

## 3. Why "the glm thing" needs a prior (penalized IRLS)

The brief was to implement these without a Stan toolchain — i.e., with base-R
GLMs. A plain `glm()` MLE works for the hurdle's two parts, but the **naive +.01
model's unpenalized MLE diverges**: with \(\sim\!45\!-\!55\%\) of the outcomes
mapped to \(\log(0.01)\), the Gamma log-link IRLS exhibits a quasi-separation —
`exp(eta)` overflows and `glm.fit` aborts with *"NA/NaN/Inf in 'x'"* (clamping
`eta` only yields degenerate runaway coefficients, e.g. \(\hat\beta_z\approx
2\times10^{13}\)). This is exactly why the paper fit these in Stan: the **prior
regularizes** the fit.

The faithful, Stan-free equivalent of a weakly-informative Bayesian GLM is a
**MAP estimate with a Laplace (large-sample Gaussian) posterior**, obtained by
**penalized IRLS** with a Gaussian prior \(N(0,\sigma_0^2)\) on each coefficient
(default \(\sigma_0=5\), weakly informative):

- point estimate: penalized IRLS update
  \(\beta\leftarrow (X'WX + \Sigma_0^{-1})^{-1} X'Wz\);
- posterior: \(\beta\sim N(\hat\beta,\ \hat\phi\,(X'WX+\Sigma_0^{-1})^{-1})\),
  with \(\hat\phi\) the Pearson dispersion (1 for the logistic part);
- posterior draws of \(\beta\) feed the G-computation, giving `cate`, `ate`,
  `p1`, `p0` (credible intervals via the across-draw quantiles used by the
  existing `calc_cate_metrics()`).

The prior plays the same role as in the paper's Stan fit; for the hurdle parts
(where the MLE already converges) it only shrinks weakly. This keeps both models
robust across all 100 seeds and all sample sizes — **no MCMC, no compilation,
and no new package dependency** (only `stats` + the already-imported
`MASS::mvrnorm`).

## 4. Outputs returned (comparability)

Both wrappers return the same objects as the other benchmarks:

- `cate` — `draws × n` matrix of posterior CATE draws,
- `ate`  — row means of `cate`,
- `p1`, `p0` — `draws × n` matrices of \(P(Y>0\mid do(a),x)\) for \(a=1,0\).

`gamma_hurdle()` returns all four (its logistic part **is** the hurdle, so the
zero-probability contrast `p1 - p0` is well defined). `gamma_plus01()` has **no
hurdle component**, so it returns `p1 = p0 = NULL`; its notebooks populate the
`*_Hurdle_*` columns with `NA` (mirroring how the BCF-Linear notebooks treat the
Smear/Hurdle columns). This matches the paper, which uses both only as
response-scale (ATE) comparators.

## 5. Package integration

- `R/gamma_benchmarks.R` — `gamma_hurdle()` and `gamma_plus01()` plus the shared
  internal `.pen_irls()` (penalized-IRLS engine), `.gamma_design()` (the
  interacted design), and `.gamma_post()` (posterior draws).
- `NAMESPACE` — `export(gamma_hurdle)`, `export(gamma_plus01)`, and
  `importFrom(stats, binomial)` / `importFrom(stats, Gamma)` (`MASS::mvrnorm`,
  `stats::plogis`, `stats::vcov` were already imported). **No** new `Imports`.

`source()`/`library(zicbcf)` then `gamma_hurdle(...)` / `gamma_plus01(...)` run
end-to-end.

## 6. Notebooks (the simulation grid)

`simulation_studies/generate_gamma_notebooks.py` generates **48** Colab
notebooks, continuing the existing `folder1..32` organisation and reusing the
**byte-identical** DGP generators and `calc_cate_metrics()` from
`generate_benchmark_notebooks.py` (so the numbers are directly comparable):

- **Gamma hurdle:** `folder33` (standard N=500) · `34/35/36` (N=1000/100/250) ·
  `37/38/39/40` (ZI levels 1/2/4/5).
- **Gamma +.01:** `folder41` (standard) · `42/43/44` (N) · `45/46/47/48` (ZI).

Each folder holds three notebooks (DGP A/B/C). Each row records a `Status`
(`OK`/`FAILED`) column and metric columns prefixed `GammaHurdle_*` /
`GammaP01_*` (mirroring the `EDP_*` / `DPglm_*` / `Linear_*` / `Smear_*`
convention). MCMC budget matches the rest of the study (100 seeds,
`NBURN = NSIM = 1000`; `NSIM` sets the number of posterior draws — `NBURN`/`NTHIN`
are accepted for interface symmetry but unused by the Laplace posterior).

## 7. Runtime

Unlike EDP (cubic per sweep, tens of hours per notebook), these GLM benchmarks
are **cheap**: a couple of penalized-IRLS fits per dataset, ~0.05–0.25 s per
model at N≤1000. A full 100-seed notebook runs in well under a minute, so the
entire Gamma grid is comfortably feasible on free Colab compute.

## 8. Validation

Local checks across DGP A/B/C and N ∈ {100, 250, 500, 1000} (8 seeds each)
returned finite metrics for every configuration. The qualitative behaviour
matches the paper: the Gamma hurdle is a reasonable parametric comparator (CATE
bias → ~0 as N grows, coverage near the nominal 0.95), while the naive Gamma
+.01 model remains biased even at N = 1000 — reproducing the paper's point that
"adding a small constant severely degrades the accuracy and precision of
treatment effect estimates."
