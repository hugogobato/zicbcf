# Zero-Inflated Continuous BCF (ZIC-BCF) Path A Implementation Plan

This plan implements **Path A (Two-Part Hurdle BCF)** as described in `ZIC_BCF_Research_Proposal.md`. 
Under the assumption of independent priors for the binary hurdle and continuous log-intensity components, the posterior factors completely:
$$p(\boldsymbol{\theta}_b, \boldsymbol{\theta}_c \mid Y, I, X, Z) = p(\boldsymbol{\theta}_b \mid I, X, Z) \cdot p(\boldsymbol{\theta}_c \mid Y^+, I = 1, X, Z)$$

To avoid naming conflicts with other paths being explored simultaneously by other agents, all files, exported C++ entry points, and R functions will be suffixed with **`_pathA`** (or **`_pathA`** casing).

---

## Proposed Naming Conventions

- **C++ Source File**: [src/zicbcf_pathA.cpp](file:///home/hugo_souto/Stuff/Research/ZI-BCF/src/zicbcf_pathA.cpp)
- **C++ Exported Function**: `zicbcfCore_pathA`
- **R Source File**: [R/zicbcf_pathA.R](file:///home/hugo_souto/Stuff/Research/ZI-BCF/R/zicbcf_pathA.R)
- **R Wrapper Function**: `zicbcf_pathA`

---

## User Review Required

> [!IMPORTANT]
> - **Unified C++ Sampler**: Although Path A's posteriors factor completely (meaning one *could* run two separate models sequentially), running them jointly in a single C++ function is highly advantageous. It guarantees that the MCMC chains are perfectly synchronized, reduces R-to-C++ calling overhead, and provides a clean, single-point API for fitting ZIC-BCF.
> - **Subpopulation Propensity Adjustment (SPA)**: For the continuous part, we will estimate a separate propensity score $\widehat{\pi}^+_i = P(Z_i=1 \mid X_i, Y_i>0)$ on the active subpopulation. Passing this to the prognostic forest of the continuous part avoids selection-induced Regularization-Induced Confounding (RIC) bias.

---

## Proposed Changes

### C++ Backend Component

#### [NEW] [src/zicbcf_pathA.cpp](file:///home/hugo_souto/Stuff/Research/ZI-BCF/src/zicbcf_pathA.cpp)
We will create a new C++ source file `src/zicbcf_pathA.cpp` that contains:
- `zicbcfCore_pathA(...)`: The exported Rcpp/Armadillo function.
- It will receive two separate data sets:
  - **Hurdle Part**: Full sample size $n$, binary indicator $y_b \in \{0, 1\}$, design matrices $X_{\mu_b}$ and $X_{\tau_b}$, moderator matrix $\Omega_{\tau_b}$, and priors.
  - **Continuous Part**: Active subset size $n_{active}$, continuous log-transformed positive outcomes $y_c$, design matrices $X_{\mu_c}$ and $X_{\tau_c}$, moderator matrix $\Omega_{\tau_c}$, and priors.
- In the MCMC loop:
  - **Hurdle (Probit BCF)**: Draw truncated latent variables $Z_i^*$ based on $y_b$ and current hurdle fit, update hurdle prognostic/moderating forests using $Z_i^*$ as a continuous outcome, and keep standard deviation $\sigma_b = 1.0$.
  - **Continuous (Gaussian BCF)**: Update continuous prognostic/moderating forests using $y_c$ as the outcome, and update continuous variance $\sigma_c^2$ using standard inverse-chi-squared draws from residuals on the active subset.
  - Save posterior fits of both components on every saving iteration.

### R Wrapper Component

#### [NEW] [R/zicbcf_pathA.R](file:///home/hugo_souto/Stuff/Research/ZI-BCF/R/zicbcf_pathA.R)
We will create a new R wrapper file `R/zicbcf_pathA.R` that provides the user-facing function `zicbcf_pathA(...)`:
- **Step 1**: Identify the active subset where $Y_i > 0$.
- **Step 2**: Estimate the subpopulation propensity score (SPA) $\widehat{\pi}^+_i = P(Z_i = 1 \mid X_i, Y_i > 0)$ on the active subset using a flexible learner (e.g. classification random forest or logistic regression).
- **Step 3**: Set up the design matrices for both the full sample (hurdle) and active subset (continuous) parts, appending the correct propensity score to the prognostic covariate matrices.
- **Step 4**: Call the C++ backend `zicbcfCore_pathA`.
- **Step 5**: Process the MCMC draws to compute potential outcomes and causal effects:
  - Hurdle probability draw: $P(Y_{di}^{(s)} > 0) = \Phi(\mu_b^{(s)} + d \cdot \tau_b^{(s)})$
  - Continuous log-mean draw: $\log Y_{di}^{+(s)} = \mu_c^{(s)} + d \cdot \tau_c^{(s)}$
  - Control mean draw: $\mu_0^{(s)}(X_i) = \Phi\left( \mu_b^{(s)} \right) \cdot \exp\left( \mu_c^{(s)} + \frac{\sigma^{2(s)}}{2} \right)$
  - Treated mean draw: $\mu_1^{(s)}(X_i) = \Phi\left( \mu_b^{(s)} + \tau_b^{(s)} \right) \cdot \exp\left( \mu_c^{(s)} + \tau_c^{(s)} + \frac{\sigma^{2(s)}}{2} \right)$
  - Response-scale CATE: $\tau_{\text{resp}}^{(s)}(X_i) = \mu_1^{(s)}(X_i) - \mu_0^{(s)}(X_i)$
- Return a structured `zicbcf_fit_pathA` list.

---

## Verification Plan

### Automated Tests
1. **Compilation**: Compile the package to ensure `zicbcf_pathA.cpp` and `zicbcf_pathA.R` compile cleanly via Rcpp/Armadillo.
2. **Simulation Study**: Set up a simulation script based on the Data-Generating Process (DGP) proposed in Section 4 of `ZIC_BCF_Research_Proposal.md`.
3. **Correctness**: Verify that:
   - Hurdle probit predictions match the true participation probabilities.
   - Continuous log-normal BCF predictions match the true log-means.
   - Response-scale CATE and ATE are accurately recovered.
