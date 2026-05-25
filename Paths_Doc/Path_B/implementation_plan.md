# Implementation Plan - Path B: Unified Single-Forest Tweedie-BCF (Isolated Version)

We propose implementing **Path B (Unified Single-Forest Tweedie-BCF)** in completely isolated files (`src/countbcf_pathb.cpp` and `R/countbcf_pathb.R`) to prevent any conflicts with other agents working on different paths.

Through a deep mathematical investigation, we have discovered that when the Tweedie power index is fixed to **$p = 1.5$** (the most standard and widely used value, corresponding to the compound Poisson-exponential process), the joint likelihood under the compound Poisson-Gamma representation is **exactly conjugate** under the existing GIG-mixture prior of `countbcf`! 

This allows us to implement the entire unified Tweedie-BCF model with **zero changes to the tree structure proposal (birth/death) or the GIG leaf parameter MCMC code**, reusing the highly optimized backfitting engine in `countbcf` by copying its structure into our new isolated files and passing transformed statistics.

---

## 1. Mathematical Framework and Conjugacy

Let the outcome $Y_i \ge 0$ follow a Tweedie distribution with mean $\mu_i$, dispersion parameter $\phi$, and power parameter $p = 1.5$:
$$Y_i \mid \mu_i, \phi \sim \text{Tweedie}(\mu_i, \phi, p = 1.5)$$

Using the compound Poisson-Gamma process representation of the Tweedie distribution:
1. $N_i \sim \text{Poisson}(\lambda_i)$ represents the latent number of Poisson events for unit $i$.
2. $X_{ij} \sim \text{Exponential}(\gamma_i)$ represents the size of event $j$ for unit $i$.
3. $Y_i = \sum_{j=1}^{N_i} X_{ij}$ is the observed continuous semicontinuous outcome.

### Parameter Mapping for $p = 1.5$
From the Tweedie analytical relationships:
- Event shape: $\alpha = \frac{2-p}{p-1} = 1.0$ (representing Exponential distribution sizes)
- Poisson rate: $\lambda_i = \frac{\mu_i^{2-p}}{\phi (2-p)} = \frac{2\sqrt{\mu_i}}{\phi}$
- Exponential scale parameter: $\gamma_i = \frac{\mu_i^{p-1}}{\phi (p-1)} = \frac{2\sqrt{\mu_i}}{\phi}$

### The Joint Likelihood & GIG Conjugacy
Conditioning on the latent event count $N_i$, the joint likelihood of $(N_i, Y_i)$ given the mean $\mu_i$ and dispersion $\phi$ is:
$$P(N_i, Y_i \mid \mu_i, \phi) = \frac{\lambda_i^{N_i} e^{-\lambda_i}}{N_i!} \frac{\gamma_i^{N_i} Y_i^{N_i - 1} e^{-\gamma_i Y_i}}{\Gamma(N_i)}$$
$$\propto (\lambda_i \gamma_i)^{N_i} \exp\left( - \lambda_i - \gamma_i Y_i \right)$$

Substituting the mappings $\lambda_i = \gamma_i = \frac{2\sqrt{\mu_i}}{\phi}$:
$$P(N_i, Y_i \mid \mu_i, \phi) \propto \left( \frac{4 \mu_i}{\phi^2} \right)^{N_i} \exp\left( - \frac{2(1 + Y_i)}{\phi} \sqrt{\mu_i} \right)$$
$$\propto \mu_i^{N_i} \exp\left( - \frac{2(1 + Y_i)}{\phi} \sqrt{\mu_i} \right)$$

Let us define the target tree function as $\theta_i = \log \sqrt{\mu_i} = 0.5 \log \mu_i$, so that $\sqrt{\mu_i} = \exp(\theta_i)$ and $\mu_i = \exp(2\theta_i)$. Substituting this into the joint likelihood yields:
$$L_i(\theta_i) \propto \exp(2 N_i \theta_i) \exp\left( - \frac{2(1 + Y_i)}{\phi} \exp(\theta_i) \right)$$

This is **mathematically identical** to the Poisson-type conjugate likelihood:
$$L_i(\theta_i) \propto \exp(u_i \theta_i) \exp\left( - v_i \exp(\theta_i) \right)$$
where:
- **Pseudo-count:** $u_i = 2 N_i$
- **Pseudo-exposure:** $v_i = \frac{2(1 + Y_i)}{\phi}$

---

## 2. Gibbs Sampling of Latent Variable $N_i$ and Dispersion $\phi$

To implement the MCMC sampler, we sample the latent counts $N_i$ and the dispersion parameter $\phi$ alongside the trees:

### Step 1: Draw Latent Variable $N_i$
- For $Y_i = 0$, $N_i = 0$ with probability $1$.
- For $Y_i > 0$, we draw $N_i \in \{1, 2, 3, \dots\}$ from its posterior:
  $$P(N_i = n \mid Y_i > 0, \mu_i, \phi) \propto \frac{ \left( \frac{4 Y_i \sqrt{\mu_i}}{\phi^2} \right)^n }{ n! (n-1)! }$$
  To ensure numerical stability and prevent overflow, we compute the term values in **log-space** iteratively:
  $$\log T_1 = \log(d_i), \quad d_i = \frac{4 Y_i \sqrt{\mu_i}}{\phi^2}$$
  $$\log T_n = \log T_{n-1} + \log d_i - \log n - \log(n-1)$$

### Step 2: Draw Dispersion Parameter $\phi$
By placing a conjugate Inverse-Gamma$(a_0, b_0)$ prior on $\phi$ (with weakly informative $a_0 = 1.0, b_0 = 1.0$), we draw $\phi$ directly using Gibbs:
$$\phi \mid N, Y, \mu \sim \text{Inverse-Gamma}\left( a_0 + 2\sum_{i: Y_i > 0} N_i, \quad b_0 + \sum_{i=1}^n 2(1 + Y_i)\sqrt{\mu_i} \right)$$

---

## 3. Proposed Changes

We will create isolated Path B files containing `countbcf_pathb` for the sampler.

### C++ Component: [NEW] [countbcf_pathb.cpp](file:///home/hugo_souto/Stuff/Research/ZI-BCF/src/countbcf_pathb.cpp)

This file will be a copy of `src/countbcf.cpp` modified as follows:
1. **Rename function:** Export `countbcf_pathb(...)` via Rcpp instead of `countbcf`.
2. **Add Helper Function:**
   Implement `draw_latent_N` at the top of the file to perform the stable log-space categorical draw of $N_i$.
3. **Initialize MCMC Variables:**
   - Define `std::vector<int> N_latent(n, 0);` to store latent counts.
   - Use the `kappa` parameter to hold the dispersion $\phi$ (so it is saved and returned seamlessly).
4. **MCMC Loop Modifications:**
   - At the start of the sweep, draw the latent counts $N_i$:
     ```cpp
     for (size_t k = 0; k < n; ++k) {
       double theta_k = log_f_per_group[0][k];
       double mu_k = exp(offset[k] + 2.0 * theta_k);
       N_latent[k] = draw_latent_N(y[k], mu_k, kappa, gen);
     }
     ```
   - Update `kappa` (dispersion $\phi$) via Gibbs:
     ```cpp
     double sum_N = 0.0;
     double sum_V = 0.0;
     for (size_t k = 0; k < n; ++k) {
       if (y[k] > 0) sum_N += N_latent[k];
       double theta_k = log_f_per_group[0][k];
       double mu_k = exp(offset[k] + 2.0 * theta_k);
       sum_V += 2.0 * (1.0 + y[k]) * sqrt(mu_k);
     }
     kappa = 1.0 / gen.gamma(1.0 + 2.0 * sum_N, 1.0 / (1.0 + sum_V));
     ```
   - Transform tree statistics in the backfitting loop:
     ```cpp
     u_vec[k]  = 2.0 * N_latent[k] * omega_k;
     r_tree[k] = (2.0 * (1.0 + y[k]) / kappa) * omega_k;
     ```
   - Ensure the tree fits update `log_mean[k]` scaled by `2.0` (since the trees predict $0.5 \log \mu_k$):
     ```cpp
     if (is_count_group) {
       log_mean[k] -= 2.0 * prior_info[s].eta * ftemp[k];
     }
     ```
     and similarly when adding it back:
     ```cpp
     if (is_count_group) {
       log_mean[k] += 2.0 * prior_info[s].eta * ftemp[k];
     }
     ```

### R Component: [NEW] [countbcf_pathb.R](file:///home/hugo_souto/Stuff/Research/ZI-BCF/R/countbcf_pathb.R)

This file will be a copy of `R/countbcf.R` modified as follows:
1. **Rename function:** Export `countbcf_pathb(...)` wrapper function.
2. **Hardcode count_model/model_type:**
   Force it to run the Tweedie code in C++ directly, using 2 forests (1 prognostic, 1 moderating), and correctly format the outputs.
3. **CATE & ATE re-transformation helper:**
   Implement response-scale ATE / CATE using the compound Poisson-exponential formulas in the return documentation and as a utility helper function.

---

## 4. Verification Plan

We will compile the package and run a simulation study mirroring the proposed Tweedie scenario in the research proposal:

1. **DGP Generation:**
   Generate $N = 1000$ points from a compound Poisson-Gamma distribution with $p=1.5$, dispersion $\phi = 1.2$, and log-mean $\log \mu_i = x_{i1} + z_i(0.5 - 0.2 x_{i2})$.
2. **Model Fit:**
   Fit `countbcf_pathb(...)`.
3. **Accuracy Checks:**
   - Verify that the Gibbs chain for the dispersion parameter $\phi$ converges close to the true value of $1.2$.
   - Confirm that the estimated response-scale potential outcomes are unbiased compared to the true conditional expectations.
