# Implementation Plan — Path C-A: Joint Copula-BCF with Path A Variance Controls

We propose **Path C-A**: a refinement of **Path E** (Collapsed Joint Copula-BCF with Active-Subset Outcome Intensity) that imports the three variance-control ingredients from **Path A** (ZIC-BCF). The motivation is the bias-variance gap observed on DGP C (Tweedie, nburn = 2000):

* Path A: CATE RMSE = **9.10** (low), |ATE Bias| = **3.67** (high)
* Path C updated: CATE RMSE = **10.98** (high), |ATE Bias| = **0.36** (low)
* Path E: CATE RMSE = **10.86**, |ATE Bias| = **1.44** — better balanced than Path C updated, but RMSE still inflated relative to Path A.

Path C-A is designed to inherit Path C updated's joint-selection ATE accuracy while pulling its CATE RMSE down to Path A territory.

---

## 1. Mathematical Diagnosis

In Path E, after data-augmenting $V_k$ for the inactive units, the joint conjugate Gibbs update for $(\beta, \sigma_0^2)$ aggregates residuals over **all** $n$ units:

$$
\sum_{i=1}^{n} \delta_i \varepsilon_i, \qquad \delta_i = W_i^{\*} - \eta_{b,i},\ \varepsilon_i = V_i - \eta_{c,i}.
$$

For inactive units, $V_k$ is sampled from $\mathcal{N}(\eta_{c,k} + \beta\,\delta_k, \sigma_0^2)$, so $\varepsilon_k = V_k - \eta_{c,k} \approx \beta\,\delta_k + \nu_k$ with $\nu_k \sim \mathcal{N}(0, \sigma_0^2)$. The contribution $\delta_k\varepsilon_k \approx \beta\,\delta_k^2 + \delta_k\nu_k$ **carries no independent information about $\beta$** — it simply reinforces the current value while injecting fresh augmentation noise. Under Tweedie-style misspecification this term amplifies CATE-level variance without buying additional bias correction.

The fix is to restrict the conjugate update to the **active** subset, where $V_i$ is observed and the residual $\varepsilon_i$ genuinely informs $(\beta, \sigma_0^2)$.

---

## 2. Three Path A Ingredients Imported into Path E

| # | Ingredient | Location | Why it matters |
|---|---|---|---|
| 1 | $(\beta, \sigma_0^2)$ Gibbs restricted to the **active** subset | `src/pathca_bcf.cpp`, Step 7 | Removes augmentation-noise amplification of leverage-1 inactive residuals. |
| 2 | Data-adaptive prior scales on the log-outcome stage: `sd_control_continuous = 2 * sd(log y[y>0])`, `sd_moderate_continuous = 0.25 * sd(log y[y>0]) / sd(z[y>0])` | `R/pathca_bcf.R` | Matches Path A's leaf-prior calibration to the empirical log-scale spread of the active subset. |
| 3 | Chipman-style $\lambda$ calibration from OLS on the active subset and SPA (active-subset propensity) for the outcome forests | `R/pathca_bcf.R` | Forecasts a realistic residual variance for the outcome stage; SPA removes selection-induced confounding inside the active subset. |

The selection forests continue to use the full-sample propensity `pihat_sel`, exactly as in Path E. Only the outcome forests use `pihat_active`.

---

## 3. File Blueprint

### C++ Core: `src/pathca_bcf.cpp`
- **Function:** `pathca_bcfCore` (entry point identical signature to `pathe_bcfCore`).
- **Key change vs Path E:** Step 7 (joint conjugate $(\beta, \sigma_0^2)$ update) iterates over `active_indices` only, using $n_c = \sum_i I_i$ instead of $n$ in both the cross-products and the posterior degrees of freedom.
- All other MCMC steps (latent $W^\*$, augmented $V$, selection forests, active outcome forests, full-sample forward prediction) remain as in Path E.

### R Wrapper: `R/pathca_bcf.R`
- **Function:** `pathca_bcf`.
- New responsibilities relative to `pathe_bcf`:
  1. Compute SPA propensity (logit GLM on active subset, predicted across all units).
  2. Build outcome covariate matrices that use `pihat_active` (selection covariates keep `pihat_sel`).
  3. Apply Path A's data-adaptive `sd_control_continuous`, `sd_moderate_continuous`, and Chipman OLS-based $\lambda$ on the active subset.

### Documentation / Namespace
- **MODIFY** `NAMESPACE`: `export(pathca_bcf)` and `export(pathca_bcfCore)`.
- **MODIFY** `R/RcppExports.R`, `src/RcppExports.cpp`: register `_countbcf_pathca_bcfCore`.

---

## 4. Test Plan

1. **Compile**: `R CMD INSTALL --library=local_lib .`
2. **Smoke test**: small synthetic Tweedie sample (N = 200, nburn = nsim = 50). Verify $\beta$ and $\sigma_0$ stay finite and posterior arrays have the right shape.
3. **Full simulation** (`simulation_studies/run_pathca.R`): DGP A, B, C with nburn = 500; DGP C also with nburn = 2000.
4. **Acceptance criteria**:
   - DGP C (nburn = 2000) CATE RMSE strictly below Path E's 10.86 and ideally close to Path A's 9.10.
   - |ATE Bias| strictly below Path A's 3.67 (ideally close to Path C updated's 0.36 / Path E's 1.44).
   - Coverage remains near nominal 95%.
