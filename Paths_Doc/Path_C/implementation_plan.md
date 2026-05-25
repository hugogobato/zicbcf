# Implementation Plan - Path C: Joint Copula-BCF with Selection on Unobservables

We propose to implement Path C of the ZIC-BCF Research Proposal: a **Joint Copula-BCF with Selection on Unobservables** (Heckman-style Selection model) using Bayesian Causal Forests.

To prevent any conflict with other agents working on different paths, all our files and core entry-points will explicitly incorporate **Path C** in their names.

## User Review Required

> [!IMPORTANT]
> - This implementation introduces a new C++ MCMC sampler (`pathc_bcfCore`) that performs joint covariance estimation of $(\beta, \sigma_0^2)$ and data augmentation of latent utilities $W_i^*$ and $V_i$.
> - We assume a conjugate Normal-Inverse-Gamma prior for the joint covariance parameters $(\beta, \sigma_0^2)$ to perform exact Gibbs updates without fragile Metropolis-Hastings sweeps.

## Proposed Changes

We will create a new C++ implementation file `src/pathc_bcf.cpp` and a corresponding R wrapper `R/pathc_bcf.R`.

### C++ Core Component

#### [NEW] [pathc_bcf.cpp](file:///home/hugo_souto/Stuff/Research/ZI-BCF/src/pathc_bcf.cpp)
We will implement the `pathc_bcfCore` function inside `src/pathc_bcf.cpp`. It will:
1. Accept four distinct BCF forests:
   - Hurdle selection prognostic forest: $\mu_b$ (`t_sel_con`)
   - Hurdle selection moderating forest: $\tau_b$ (`t_sel_mod`)
   - Continuous intensity outcome prognostic forest: $\mu_c$ (`t_out_con`)
   - Continuous intensity outcome moderating forest: $\tau_c$ (`t_out_mod`)
2. Maintain latent vectors:
   - $W_i^*$ (latent utility for selection)
   - $V_i$ (augmented log-intensity, where $V_i = \log(Y_i)$ for $Y_i > 0$ and is augmented for $Y_i = 0$)
3. Run MCMC loops:
   - **Step 1:** Sample $W_i^*$ using truncated normal draws:
     - For $Y_i > 0$: $W_i^* \sim \mathcal{N}\left( \eta_{b, i} + \frac{\beta}{\sigma^2}(V_i - \eta_{c, i}), 1 - \rho^2 \right)$ truncated to $W_i^* > 0$.
     - For $Y_i = 0$: $W_i^* \sim \mathcal{N}\left( \eta_{b, i}, 1 \right)$ truncated to $W_i^* \le 0$.
   - **Step 2:** Sample $V_i$ for censored units ($Y_i = 0$):
     - $V_i \sim \mathcal{N}\left( \eta_{c, i} + \beta(W_i^* - \eta_{b, i}), \sigma_0^2 \right)$.
   - **Step 3:** Update Selection forests ($\mu_b, \tau_b$) using standard linear BCF updates with:
     - Residual target response: $Y^{(b)}_i = W_i^* - \frac{\beta}{\sigma^2}(V_i - \eta_{c, i})$.
     - Effective residual standard deviation: $\sigma_{sel} = \sqrt{1 - \rho^2}$.
   - **Step 4:** Update Outcome forests ($\mu_c, \tau_c$) using standard linear BCF updates with:
     - Residual target response: $Y^{(c)}_i = V_i - \beta(W_i^* - \eta_{b, i})$.
     - Effective residual standard deviation: $\sigma_{out} = \sigma_0 = \sqrt{\sigma_0^2}$.
   - **Step 5:** Sample covariance parameters $(\beta, \sigma_0^2)$ jointly using standard Normal-Inverse-Gamma conjugate updates:
     - Residuals: $\delta_i = W_i^* - \eta_{b, i}$, $\epsilon_i = V_i - \eta_{c, i}$.
     - Regression: $\epsilon_i = \beta \delta_i + \eta_i, \quad \eta_i \sim \mathcal{N}(0, \sigma_0^2)$.
     - Draw $\sigma_0^2$ from $\text{Inv-Gamma}(a_N, b_N)$.
     - Draw $\beta$ from $\mathcal{N}(M_\beta, \sigma_0^2 V_\beta)$.
     - Compute $\sigma^2 = \sigma_0^2 + \beta^2$ and $\rho = \beta / \sqrt{\sigma_0^2 + \beta^2}$.
   - **Step 6:** Update the PX parameters `eta` for each of the 4 forests.

---

### R Interface Component

#### [NEW] [pathc_bcf.R](file:///home/hugo_souto/Stuff/Research/ZI-BCF/R/pathc_bcf.R)
We will write an R interface function `pathc_bcf` that:
1. Performs input verification and validation on covariates and outcomes.
2. Standardizes outcomes and computes initial estimates of the propensity scores.
3. Automatically formats treatment/covariate matrices for the C++ backend.
4. Invokes the generated `pathc_bcfCore` function.
5. Returns a structured list of MCMC draws for both hurdle and intensity parameters, including $\rho$ and $\sigma^2$.

---

## Verification Plan

### Automated Tests
We will write a test script `/home/hugo_souto/Stuff/Research/ZI-BCF/example/run_pathc_bcf.R` to:
1. Simulate data from the Heckman selection data generating process (DGP).
2. Fit our new `pathc_bcf` model.
3. Verify that the MCMC chain runs without crashes and successfully converges to estimate the selection correlation $\rho$ and treatment effects on both selection and outcome scales.
