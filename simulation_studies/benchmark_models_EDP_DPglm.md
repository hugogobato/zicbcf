# Integration of the EDP+ZI and DPglm Benchmark Models

*Documentation of how the two Bayesian nonparametric benchmark models were
obtained, adapted, and integrated into the `zicbcf` package and the simulation
study. Written to support the methods / reproducibility sections of the paper.*

Last updated: 2026-06-01.

---

## 1. Objective

The first phase of the simulation study compared two proposed estimators —
**BCF-Linear** (`bcf_continuous_linear`) and **ZIC-BCF-Smear** (`zicbcf_smear`) —
across data-generating processes (DGPs) A, B, C and a battery of sample-size and
zero-inflation sensitivity settings.

To benchmark these against the existing Bayesian nonparametric (BNP) literature
for zero-inflated causal effects, we add the two competing methods used by
**Kim, Li, Xu & Liao (2024)**, *"Bayesian nonparametric model for heterogeneous
treatment effects with zero-inflated data"*, *Statistics in Medicine*
43(30):5968–5982 (doi:10.1002/sim.10266):

1. **EDP+ZI** — the enriched Dirichlet process zero-inflated mixture (the
   *proposed* method of that paper).
2. **DPglm** — a Dirichlet-process mixture of zero-inflated regressions
   (Oganisian et al., 2021), used there as a competing method.

Goal: run both benchmarks on the **same** DGPs / sensitivity grid as the
proposed estimators, with results directly comparable, and make them runnable on
Google Colab through a single `devtools::install_github("hugogobato/zicbcf")`.

## 2. Source of the reference code

The reference implementation was taken from the authors' public repository
`https://github.com/lit777/EDPcausal` (the local copy was the `EDPcausal-main/`
folder, scenarios `sim1/` and `sim2/`). The relevant files:

| Model | Files used | Language |
|-------|-----------|----------|
| EDP+ZI | `sim2/clustering.cpp`, `sim2/MCMC_utils.h`, `sim2/function.R`, driver `sim2/sim_edp.R` | RcppArmadillo + R |
| DPglm  | `sim1/SourceCode/DPMix_Full_SourceCode_cpp_new.R`, `sim1/SourceCode/class_update.R`, driver `sim1/sim_dpglm.R` | R |

`sim1/clustering.cpp` and `sim1/function.R` are byte-identical to the `sim2/`
versions. The DPglm `GcompSupport.cpp` / `G_comp_Vectorized_*.R` files are **not**
needed: the driver uses `DPglmMix`'s in-sample posterior-predictive draws
(`pp$y_isp1/0`, `pp$z_isp1/0`), not the separate G-computation helpers.

After integration the `EDPcausal-main/` folder is no longer referenced by
anything and can be deleted.

## 3. The two models (what they compute)

Both target the conditional average treatment effect (CATE)
\(\tau(x)=E[Y(1)-Y(0)\mid X=x]\) and the zero-probability contrast
\(\tau_\pi(x)=P(Y=0\mid do(1),x)-P(Y=0\mid do(0),x)\) for a **zero-inflated /
hurdle** outcome \(Y\ge 0\): a point mass at 0 mixed with a positive continuous
part. Both place a (nested) Dirichlet-process prior over latent clusters; within
a cluster the outcome is a zero-inflated truncated-normal with
\(P(Y=0)=\mathrm{expit}(r)\) and the positive part \(\sim N_{(0,\infty)}(x\beta,\sigma^2)\).
EDP+ZI uses an **enriched** DP (covariate clusters nested within outcome
clusters); DPglm uses a single DP over the joint \((Y,A,X)\).

For comparability with the proposed estimators, both wrappers return the same
objects:

- `cate` — `draws × n` matrix of posterior CATE draws,
- `ate`  — row means of `cate`,
- `p1`, `p0` — `draws × n` matrices of \(P(Y>0\mid do(a),x)\) for \(a=1,0\)
  (so the hurdle contrast is `p1 - p0`).

These feed the **identical** `calc_cate_metrics()` used by the existing
notebooks (RMSE, |bias|, 95% coverage, correlation, CI length, estimated ATE).

## 4. Why an adaptation was unavoidable (covariate structure)

The authors' simulation generates **one binary covariate \(X_1\), one discrete
\(X_2\), and three continuous \(X_3\!-\!X_5\)**, and their EDP C++ code is
*hard-wired* to that layout: the design matrix is
`xa = (1, A, X1[binary], X2, X3, X4, X5)`, with `A` and `X1` modelled by
categorical (Dirichlet/Bernoulli) kernels and `X2–X5` by Gaussian kernels.

Our DGPs A/B/C generate **five continuous covariates** \(X_1\!-\!X_5\sim N(0,1)\)
plus a binary treatment \(Z\). Feeding a continuous \(X_1\) into the binary
kernel is invalid (the `dbinom` likelihood is undefined off \(\{0,1\}\)).

**Adaptation (EDP only):** the treatment keeps its categorical kernel and **all
five covariates use the Gaussian kernel** (`xPAR_mu`/`xPAR_sig` extended from 4
to 5 rows; the binary-\(X_1\) kernel `xPAR_p1` removed). This is the model's own
prescribed handling of continuous covariates — the paper states covariates "can
be continuous or categorical" and already applies Gaussian kernels to the
continuous ones — so it changes the *configuration*, not the model.

**DPglm needs no such change:** `DPglmMix` is parameterized by `x_type`, so we
simply pass `x_type = c("numeric"×5, "binary")`.

## 5. Implementation defects found in the reference code, and how they were handled

The published **model** is sound and peer-reviewed; the published **code** is
research-grade and contains genuine defects. Per the project rule — *implement
the paper's model exactly, but fix implementation bugs and numerical defects* —
the model (kernels, priors, links, DP weights, G-computation) is untouched, and
the following code-level defects are corrected. Each correction was verified to
change only the numerics, never the model.

### EDP+ZI (`clustering.cpp` / `function.R`)

1. **`1/e` constant instead of the expit (genuine bug).** The cluster-assignment
   likelihood wrote the structural-zero probability as `exp(r)/exp(1+r)`, which
   algebraically equals \(e^{r-1-r}=e^{-1}\approx0.368\) — a constant independent
   of covariates and parameters. The model specifies \(\pi=\mathrm{expit}(r)=
   \exp(r)/(1+\exp(r))\) (and the parameter updates and post-processing in the
   same codebase *do* use the correct expit). Restored to the model's expit,
   evaluated stably as \(-\log1p(e^{\mp r})\).

2. **Frozen concentration-parameter update.** `update_parameters` read
   `alpha_theta`/`alpha_omega` from the global scope, where they were fixed at
   their initial value `2` and never advanced, so the Escobar–West augmentation
   restarted from `2` every iteration. Corrected to thread the *current* values
   under the paper's Gamma(1,1) priors (a proper update; the prior and
   augmentation scheme are unchanged).

3. **Product-space under/overflow → fatal aborts.** Cluster membership
   probabilities were accumulated as a product of Gaussian/zero-inflated kernels.
   For an outlier subject the product underflows to exactly 0, making the whole
   weight vector zero and Rcpp's `sample()` throw *"Too few positive
   probabilities"*. The authors' driver hides this by wrapping every replicate in
   `tryCatch(..., error=function(e){})`, silently **discarding** such datasets.
   We instead accumulate the **same** weights in log-space and exp-normalise
   (subtracting the max), which is the standard, underflow-free implementation of
   the identical Algorithm-8 weights.

4. **Truncated-normal kernel via an R callback.** The outcome kernel was
   evaluated by calling `truncnorm::dtruncnorm()` from C++ once per cluster per
   subject (slow R callbacks inside the inner loop). Replaced by the closed-form
   log density \(\log N(y\mid x\beta,\sigma)-\log\Phi(x\beta/\sigma)\) on
   \((0,\infty)\) — exact, and removes the callback.

5. **Post-processing NaNs (numerical).** The G-computation membership weight
   `prob <- prob/rowSums(prob)` produces `0/0 = NaN` for outliers (all kernels
   underflow); and the per-cluster expectation \(E[Y_{>0}]\) via
   `truncnorm::etruncnorm` returns `NaN`/`Inf` for near-degenerate clusters (tiny
   \(\sigma\): the Mills ratio \(\phi(a)/(1-\Phi(a))\) suffers catastrophic
   cancellation). Both are computed as their exact mathematical limits:
   the membership via log-sum-exp softmax, and the truncated mean via the same
   closed form with exact clamping of the extreme-`a` limits (\(a\gg0\to0\),
   \(a\ll0\to m\)). The authors' Monte-Carlo `mean(rtruncnorm(1000,…))` is
   robust to the tiny-\(\sigma\) case but is a noisy, far slower estimator of the
   exact value we compute.

The original yPAR_sig "vector indexed as a matrix" quirk (`yPAR_sig[j,1]`) — which
only worked because `wrap(arma::vec)` returns an \(n\times1\) matrix — is made
explicit in the rewrite.

### DPglm (`DPMix_Full_SourceCode_cpp_new.R` / `class_update.R`)

The DPglm sampler is cleaner and was kept **verbatim** (it is the authors'
exact code, included unmodified as internal functions). Only a thin wrapper
`dpglm()` was added to: build the design (`x_type` all-numeric + binary
treatment), supply the authors' OLS-based \(\beta\) prior and their
`Normal(0,2)` / `Inv-Gamma(10,10000)` hyper-priors, and reshape the in-sample
posterior-predictive draws into `cate`/`ate`/`p1`/`p0`. No model change.

## 6. Numerical-equivalence checks

For datasets where the authors' code runs without aborting and without producing
NaNs, the corrected computations agree with the reference quantities to floating
point (e.g. the corrected truncated-mean equals `etruncnorm` wherever the latter
is finite). The corrections only *add* well-defined values where the reference
code aborted or returned NaN/Inf. Local validation (DGP A, several seeds) shows
finite output everywhere and the expected qualitative behaviour reported in the
paper (DPglm shows clear ATE bias; EDP+ZI is less biased but variable).

## 7. Package integration

Both benchmarks now ship **inside the `zicbcf` package**:

- `src/edp_clustering.cpp` — adapted/corrected EDP cluster sampler, exported via
  Rcpp as `edp_clustering_c` (registered through `Rcpp::compileAttributes()`).
- `R/edp_zi.R` — `edp_zi()`: initial state, MCMC loop, the (adapted)
  `update_parameters`, and the log-stable G-computation.
- `R/dpglm.R` — `dpglm()` wrapper plus the authors' `DPglmMix` / `class_update`
  source (verbatim, internal).
- `DESCRIPTION` — added `Imports: stats, truncnorm, mvtnorm, MASS,
  LaplacesDemon, invgamma`.
- `NAMESPACE` — `export(edp_zi)`, `export(dpglm)`, plus the required
  `importFrom` lines.

`R CMD INSTALL` succeeds and `library(zicbcf); edp_zi(...); dpglm(...)` run
end-to-end from the installed package.

## 8. Notebooks (the simulation grid)

`simulation_studies/generate_benchmark_notebooks.py` generates **48** Colab
notebooks (mirroring the existing `colab_notebooks/folder1–16` organisation):

- **EDP+ZI:** `folder17` (standard N=500) · `18/19/20` (N=1000/100/250) ·
  `21/22/23/24` (ZI levels 1/2/4/5).
- **DPglm:** `folder25` (standard) · `26/27/28` (N) · `29/30/31/32` (ZI).

Each folder holds three notebooks (DGP A/B/C). The DGP generators and
`calc_cate_metrics()` are **byte-identical** to the existing notebooks, so the
benchmark numbers are directly comparable. MCMC budget matches the existing
study: 100 seeds, `nburn = nsim = 1000`. Each row records a `Status`
(`OK`/`FAILED`) column and metric columns prefixed `EDP_*` / `DPglm_*`
(mirroring the existing `Linear_*` / `Smear_*` convention), so a failed seed is
recorded transparently rather than silently dropped.

## 9. Runtime caveat (important)

The EDP sampler is inherent **O(n³) per MCMC sweep** (a pairwise duplicate-row
scan of the cluster labels is redone for every subject, plus a per-subject R
`order()` callback) — this is a property of the authors' algorithm/implementation,
not of the adaptation. Indicative local timing: ~0.16 s/sweep at n=250, so a full
EDP notebook at **N=500** (2000 sweeps × 100 seeds) is on the order of tens of
hours, and the **N=1000** notebooks are several times worse. DPglm (pure-R Gibbs)
is also far slower than the BCF estimators but markedly cheaper than EDP.
Practical options if this is prohibitive on Colab: reduce the number of seeds for
EDP, skip the N=1000 EDP cells, or optimise the O(n³) duplicate-row scan to
O(n log n) (a model-preserving speed-up that can be added on request).

## 10. Reproducibility summary

1. `devtools::install_github("hugogobato/zicbcf")` (pulls in all benchmark deps).
2. Open any `colab_notebooks/folder17–32` notebook and run all cells.
3. Each notebook writes `results_{edp|dpglm}_{dgp}_{N}.csv`.
4. `EDPcausal-main/` is no longer needed and can be removed.
