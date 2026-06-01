# Implementation Plan - Path E: Collapsed Joint Copula-BCF (Active Outcome Fit)

We propose **Path E**: a **Collapsed Joint Copula-BCF with Active-Subset Outcome Intensity**.

Path E resolves the CATE RMSE inflation of joint selection models by grow/pruning and leaf-updating the outcome forests only on the active subset ($Y_i > 0$), while still updating the selection forests and joint covariance parameters $(\beta, \sigma_0^2)$ on the full sample.

---

## 1. Mathematical Design

By integrating out the unobserved outcomes $V_i$ for inactive units ($I_i = 0$), the outcome forest trees and leaf parameter updates can be collapsed to the active subset.

For active units, the target is adjusted for selection bias:
$$\text{Target}_i = V_i - \beta(W_i^* - \eta_{b,i})$$
The outcome forests are updated only on these units using standard residual variance $\sigma_0^2$. At the end of each sweep, we perform a fast forward prediction pass using the full design matrices to obtain $\eta_{c, i}$ for all units (required for data-augmenting selection utilities and updating the joint covariance).

---

## 2. File Blueprint

### C++ Core: `src/pathe_bcf.cpp`
- **Function:** `pathe_bcfCore`
- **Data Wrappers:**
  - `di_out_con` and `di_out_mod` of size $n_{\text{active}}$, used for forest grow/prune and leaf updates.
  - `di_out_con_full` and `di_out_mod_full` of size $n$, used for fast forward prediction at the end of each iteration.
- **MCMC Loop:**
  1. Draw latent selection utilities $W_i^*$ for all $n$ units.
  2. Draw augmented log-intensity $V_i$ for inactive units ($I_i = 0$).
  3. Update selection forests ($\mu_b, \tau_b$) on the full sample.
  4. Update outcome forests ($\mu_c, \tau_c$) on the active sample only.
  5. Predict outcome forest fits $\eta_{c, i}$ for all $n$ units.
  6. Joint conjugate Gibbs update for covariance $(\beta, \sigma_0^2)$.
  7. Save unmodulated treatment effects $\tau_b(X_i)$ and $\tau_c(X_i)$ for all units.

### R Wrapper: `R/pathe_bcf.R`
- **Function:** `pathe_bcf`
- **Logic:**
  1. Partitions outcome covariates and propensity scores into full and active subsets.
  2. Constructs basis matrices for both active outcome forests and full outcome prognostic/moderating designs.
  3. Calls `pathe_bcfCore` and scales resulting draws back to the log-outcome scale.
