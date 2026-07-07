# zicbcf

**Zero-Inflated Continuous Bayesian Causal Forests (ZIC-BCF) for Semicontinuous Causal Inference**

The `zicbcf` package implements state-of-the-art Bayesian Causal Forest (BCF) models optimized for semicontinuous (zero-inflated continuous) and standard continuous/binary outcomes. The core contribution is **ZIC-BCF-Smear** (ZIC-BCF with Duan's Smearing Re-transformation), which resolves scale-misspecification biases under heavily skewed or heteroskedastic semicontinuous processes (e.g. Tweedie or Gamma processes).

---

## 1. Supported Models

1. **`bcf_continuous_linear`**: Standard continuous BCF fitted directly on continuous outcomes.
2. **`bcf_binary`**: Standard binary BCF for binary outcomes.
3. **`zicbcf_smear`**: Two-part Hurdle BCF paired with Duan's Non-Parametric Smearing Re-transformation.
   * **Hurdle Stage**: Modeled via a Probit BCF on the full sample.
   * **Continuous Stage**: Modeled via a log-Gaussian BCF strictly on the active subset ($Y > 0$) with Subpopulation Propensity Adjustment (SPA) to protect against subpopulation confounding.
   * **Re-transformation**: Uses Duan's smearing factor calculated from empirical log-scale active residuals to re-transform potential outcomes back to the raw response scale, avoiding parametric misspecification bias.

---

## 2. Installation

The package compiles a high-performance C++ backend (Rcpp, RcppArmadillo, Rcereal). To install the package in R:

```R
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("hugogobato/zicbcf")
library(zicbcf)
```

### Dependencies
The package automatically installs:
* `Rcpp`
* `RcppArmadillo`
* `Rcereal`

---

## 3. Quick Usage Example

Below is a minimal example showing how to fit the proposed `zicbcf_smear` model and extract CATE and ATE draws:

```R
library(zicbcf)

# 1. Generate semicontinuous data (e.g., Log-Normal Hurdle)
set.seed(123)
n <- 500
p <- 5
X <- matrix(rnorm(n * p), n, p)

# Propensity score and treatment assignment
pihat <- plogis(-0.5 + 0.3 * X[, 1] + 0.2 * X[, 2])
z <- rbinom(n, 1, pihat)

# Hurdle indicators and continuous intensity
p_hurdle <- pnorm(0.2 + 0.5 * X[, 1] - 0.3 * X[, 3] + z * (0.4 + 0.2 * X[, 1]))
I <- rbinom(n, 1, p_hurdle)
y_pos <- exp(1.5 + 0.8 * X[, 2] + z * (0.5 - 0.3 * X[, 2]) + rnorm(n, 0, 0.5))
y <- I * y_pos

# 2. Fit ZIC-BCF-Smear
fit <- zicbcf_smear(
  y = y,
  z = z,
  x_control = X,
  x_moderate = X,
  pihat = pihat,
  nburn = 500,
  nsim = 1000
)

# 3. Extract CATE and ATE draws
cate_draws <- fit$cate  # nsim x n matrix
ate_draws <- fit$ate    # nsim vector of ATE draws

# Summarize results
cate_estimates <- colMeans(cate_draws)  # Posterior mean CATE per unit
ate_estimate <- mean(ate_draws)         # Posterior mean ATE
ate_ci <- quantile(ate_draws, c(0.025, 0.975)) # 95% Credible Interval
```

### Hurdle (Participation) Treatment Effect

Because `zicbcf_smear` is a two-part model, it also estimates the **hurdle-margin treatment effect**: the causal effect of treatment on the *probability of a positive (non-zero) outcome*, separately from the overall dollar-scale CATE. The probit hurdle stage returns the prognostic (`fit$mu_b`) and treatment-effect (`fit$tau_b`) draws on the latent probit scale (each an `nsim x n` matrix). Map them to control/treatment participation probabilities and difference them to obtain the hurdle CATE:

```R
# Participation probabilities under control (p0) and treatment (p1)
p0_draws <- pnorm(fit$mu_b)               # nsim x n : Pr(y > 0 | z = 0)
p1_draws <- pnorm(fit$mu_b + fit$tau_b)   # nsim x n : Pr(y > 0 | z = 1)

# Hurdle-margin CATE = effect on the probability of any positive outcome
hurdle_cate_draws <- p1_draws - p0_draws  # nsim x n matrix
hurdle_ate_draws  <- rowMeans(hurdle_cate_draws)  # nsim vector of hurdle ATE draws

# Summarize the hurdle (participation) effect
hurdle_cate_estimates <- colMeans(hurdle_cate_draws)   # Posterior mean hurdle CATE per unit
hurdle_ate_estimate   <- mean(hurdle_ate_draws)        # Posterior mean hurdle ATE (probability points)
hurdle_ate_ci         <- quantile(hurdle_ate_draws, c(0.025, 0.975)) # 95% Credible Interval
```

Subgroup hurdle CATEs and between-group contrasts are obtained by averaging `hurdle_cate_draws` over the units in each subgroup *within every posterior draw*, exactly as for the dollar-scale CATE. See `applied_study/run_zicbcf_hurdle_subgroups.R` for a complete worked example on the MEPS dental-expenditure data.

---

## 4. Package Internals and C++ Engine

* The core Gibbs/Metropolis-within-Gibbs samplers are implemented in C++ for maximum execution speed.
* Rcpp modules serialize the BART tree structures, enabling post-processing diagnostics via `tree_samples`.
* Subpopulation Propensity Adjustment (SPA) estimates the selection propensity dynamically on the active subset using logistic regression.

---

## 5. Simulation Studies

The `simulation_studies/` directory contains replication code and results comparing `zicbcf` to standard continuous models and joint selection models across different Data Generating Processes (DGPs).
> [!NOTE]
> The simulation studies and paper are currently a work in progress.

---

## Author

**Hugo Gobato Souto**  
Dell Technologies  
[hugo.souto@dell.com](mailto:hugo.souto@dell.com)  
ORCID: [0000-0002-7039-0572](https://orcid.org/0000-0002-7039-0572)

---

## License

GPL (>= 3)
