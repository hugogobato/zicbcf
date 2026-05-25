# Implementation Plan - Path D: Gamma Hurdle BCF via GIG Conjugacy (Path-Specific Naming)

This implementation plan outlines the creation of the C++ and R code to support **Path D: Gamma Hurdle BCF via GIG Conjugacy** as proposed in the Research Proposal.

Path D uses a two-part hurdle model where:
1. The participation (hurdle) stage is modeled via a Probit BCF on the full sample. This is already supported by the package using existing probit BCF structures.
2. The continuous intensity stage is modeled via a Gamma distribution with a log-link and GIG-conjugate leaf updates, fitted **only on the active (positive) subpopulation** ($Y_i > 0$), utilizing the **Subpopulation Propensity Adjustment (SPA)**.

To avoid collisions with other agents working on different paths concurrently, all new code files are strictly named with the suffix/prefix **`pathd_`**.

---

## Proposed Changes

We will create a new C++ file `src/pathd_gammabcf.cpp` for the C++ backend and a new R file `R/pathd_gammabcf.R` for the R wrapper.

### 1. `src` Component

#### [NEW] [pathd_gammabcf.cpp](file:///home/hugo_souto/Stuff/Research/ZI-BCF/src/pathd_gammabcf.cpp)
This file will contain the Rcpp-exported C++ function `pathd_gammabcf` which runs the MCMC sampler for the log-linear Gamma BCF model.
- **Likelihood & GIG Equivalence**:
  For active units, the outcomes are modeled as:
  $$Y_i^+ \sim \text{Gamma}\left(\kappa_c, \frac{\lambda_i}{\kappa_c}\right) \implies \log \lambda_i = \mu_c(X_i) + Z_i \tau_c(X_i)$$
  Letting $\theta_i = -\log \lambda_i$, we obtain a Poisson-type conjugate likelihood:
  $$L_i \propto \exp(\theta_i)^{\kappa_c} \exp\left( - (\kappa_c y_i) \exp(\theta_i) \right)$$
  which is mathematically identical to a Poisson model with pseudo-count $u_i = \kappa_c$ and pseudo-exposure $v_i = \kappa_c y_i$.
- **MCMC Steps**:
  - Sample the global Gamma shape parameter $\kappa_c$ using a random-walk Metropolis-Hastings step with a beta-prime prior on $\frac{\kappa_c}{1 + \kappa_c}$ (matching the prior structure for the NB dispersion parameter in `countbcf`).
  - Update prognostic ($\mu_c$) and moderating ($\tau_c$) trees using GIG-conjugate backfitting:
    - Subtract the current tree's fit.
    - Set the sufficient statistics:
      - $u_i = \kappa_c \omega_i$
      - $v_i = \kappa_c y_i \exp(\theta_{-j, i}) \omega_i$
      where $\omega_i = 1$ for prognostic trees and $\omega_i = Z_i$ for moderating trees.
    - Propose tree structural changes using `bd_loglinear` and sample leaf values using `drmu_loglinear`.
    - Add the new tree's fit back.

### 2. `R` Component

#### [NEW] [pathd_gammabcf.R](file:///home/hugo_souto/Stuff/Research/ZI-BCF/R/pathd_gammabcf.R)
An R wrapper `pathd_gammabcf` that interfaces with the compiled C++ function.
- Performs input validation and checks.
- Sets up default hyperparameter values for the GIG mixture leaf prior (concentration parameters `a0` and `a0_tau`).
- Constructs the list of specs and designs for the prognostic ($\mu_c$) and moderating ($\tau_c$) forests.
- Invokes `.Call("_countbcf_pathd_gammabcf", ...)` to run the C++ MCMC sampler.
- Returns posterior draws of the log-mean, per-forest fits, per-forest raw coefficients, and the shape parameter $\kappa_c$.

---

## Verification Plan

### Automated Compilation and Tests
Since Rcpp relies on `compileAttributes()`, the user will compile and install the package:
```R
Rcpp::compileAttributes()
devtools::document()
devtools::install()
```

### Simulation Study Verification
We can write a test simulation script in the workspace to verify that:
1. The new C++ code compiles successfully without warnings.
2. The `pathd_gammabcf` function successfully runs MCMC sweeps on simulated Gamma data.
3. The shape parameter $\kappa_c$ converges to the true data-generating value.
4. CATE estimates on the response scale show correct calibration and low RMSE.
