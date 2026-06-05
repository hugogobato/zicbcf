# Updated Semicontinuous Causal Simulation Study Framework

This report documents the updated framework for the simulation studies comparing **Standard Gaussian BCF (fitted directly on raw Y)** and **ZIC-BCF-Smear** (a two-part hurdle model with Duan's Non-Parametric Smearing Re-transformation). 

Following your updated instructions, the simulation framework has been adjusted as follows:
1. **Model Selection**: Focused exclusively on Standard Continuous BCF (Linear) and ZIC-BCF-Smear.
2. **DGP B Replacement**: Replaced the previous Gaussian (symmetric) Hurdle DGP B with a **Gamma Hurdle DGP** to represent zero-inflated and highly right-skewed semicontinuous outcomes that are distinct from Log-Normal (DGP A) and Tweedie (DGP C).
3. **Treatment Effect Scaling**: Halved the treatment effect level ($0.5$ multiplier) for DGP A and DGP C, and calibrated the new Gamma Hurdle DGP B to have a similar treatment effect magnitude.
4. **Standard Simulation Settings**: 
   - Standard sample size $N = 500$ (instead of $N = 1000$).
   - Number of simulations (seeds) $N_{\text{sim}} = 100$.
   - MCMC parameters: $N_{\text{burn}} = 1000$ and $N_{\text{sim}} = 1000$ (optimized for convergence and Colab efficiency).
5. **Colab-Ready Notebooks**: Created 48 modular self-contained Jupyter notebooks with R kernels (one per model, per level, per DGP) that install `zicbcf` from GitHub, execute the simulations, and output a single clean CSV containing the exact metrics specified in `run_5seed_single.R`.

---

## 1. Updated Semicontinuous Data-Generating Processes (DGPs)

For all three DGPs, the simulation generates $P = 5$ baseline covariates $X \sim \mathcal{N}(0, I_5)$ and a confounded treatment propensity $\pi_i = \Phi(-0.5 + 0.4 X_{i1} + 0.3 X_{i2}^2)$ with binary assignment $Z_i \sim \text{Bernoulli}(\pi_i)$.

### DGP A: Log-Normal Hurdle DGP (Skewed, Zero-Inflated)
* **Binary Hurdle:** Probit probability of positive outcome:
  $$P(I_i = 1 \mid X_i, Z_i) = \Phi\left( c_i + 0.5 X_{i1} - 0.3 X_{i3} + Z_i (0.2 + 0.1 X_{i1}) \right)$$
  *(Standard hurdle intercept $c_i = 0.2$ yields $\sim 37\%$ zeros).*
* **Continuous Intensity:** Lognormal positive outcome:
  $$\log(Y_i^+) \sim \mathcal{N}\left( 1.5 + 0.8 X_{i2} + 0.4 X_{i4} + Z_i (0.25 - 0.15 X_{i2}), 0.25 \right)$$
* **Observed Outcome:** $Y_i = I_i \cdot Y_i^+$.
* **True CATE (Response Scale):**
  $$\tau(X_i) = \Phi(c_i + 0.5 X_{i1} - 0.3 X_{i3} + 0.2 + 0.1 X_{i1}) \cdot e^{1.625 + 0.8 X_{i2} + 0.4 X_{i4} - 0.15 X_{i2}} - \Phi(c_i + 0.5 X_{i1} - 0.3 X_{i3}) \cdot e^{1.625 + 0.8 X_{i2} + 0.4 X_{i4}}$$

### DGP B: Gamma Hurdle DGP (New Replacement, Skewed, Zero-Inflated)
* **Binary Hurdle:** Same probit structure as DGP A.
* **Continuous Intensity:** Gamma positive outcome:
  $$Y_i^+ \sim \text{Gamma}\left(\text{shape} = \alpha, \text{scale} = \mu_{c, i} / \alpha \right) \quad \text{with } \alpha = 2.0$$
  where the active log-mean is designed to mirror DGP A:
  $$\log(\mu_{c, i}) = 1.5 + 0.8 X_{i2} + 0.4 X_{i4} + Z_i (0.25 - 0.15 X_{i2})$$
* **Observed Outcome:** $Y_i = I_i \cdot Y_i^+$.
* **True CATE (Response Scale):**
  $$\tau(X_i) = \Phi(c_i + 0.5 X_{i1} - 0.3 X_{i3} + 0.2 + 0.1 X_{i1}) \cdot e^{1.75 + 0.8 X_{i2} + 0.4 X_{i4} - 0.15 X_{i2}} - \Phi(c_i + 0.5 X_{i1} - 0.3 X_{i3}) \cdot e^{1.5 + 0.8 X_{i2} + 0.4 X_{i4}}$$

### DGP C: Tweedie Compound Poisson-Gamma DGP (Intrinsic Semicontinuous)
* **True Log-Mean Parameter:**
  $$\log \mu_i = 1.2 + c_i + 0.8 X_{i1} - 0.4 X_{i3} + Z_i (0.3 + 0.15 X_{i1})$$
  *(Standard mean shift $c_i = 0.0$ yields $\sim 50\%$ zeros).*
* **Outcome Generation:** Compound Poisson-Gamma Tweedie process:
  $$Y_i \sim \text{Tweedie}\left(\mu_i, \phi=1.5, p=1.5\right)$$
* **True CATE (Response Scale):**
  $$\tau(X_i) = \exp\left( 1.2 + c_i + 0.8 X_{i1} - 0.4 X_{i3} + 0.3 + 0.15 X_{i1} \right) - \exp\left( 1.2 + c_i + 0.8 X_{i1} - 0.4 X_{i3} \right)$$

---

## 2. Organization of Colab Jupyter Notebooks

All Jupyter notebooks are located in the [colab_notebooks](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks) directory. They are designed to run a single model, level, and DGP combination at a time, outputting only a CSV containing the simulation results to optimize Google Colab's execution speeds.

### 2.1 Standard Main Runs ($N = 500$, Standard Zero-Inflation)
1. **[run_bcf_dgp_a.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_bcf_dgp_a.ipynb)**: Standard BCF on DGP A
2. **[run_bcf_dgp_b.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_bcf_dgp_b.ipynb)**: Standard BCF on DGP B
3. **[run_bcf_dgp_c.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_bcf_dgp_c.ipynb)**: Standard BCF on DGP C
4. **[run_zicbcf_smear_dgp_a.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_zicbcf_smear_dgp_a.ipynb)**: ZIC-BCF-Smear on DGP A
5. **[run_zicbcf_smear_dgp_b.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_zicbcf_smear_dgp_b.ipynb)**: ZIC-BCF-Smear on DGP B
6. **[run_zicbcf_smear_dgp_c.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_zicbcf_smear_dgp_c.ipynb)**: ZIC-BCF-Smear on DGP C

### 2.2 Sample Size (N) Sensitivity Analysis ($N \in \{100, 250, 1000\}$)
- **Standard BCF (Linear)**:
  7. **[run_bcf_sensitivity_n100_dgp_a.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_bcf_sensitivity_n100_dgp_a.ipynb)**
  8. **[run_bcf_sensitivity_n100_dgp_b.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_bcf_sensitivity_n100_dgp_b.ipynb)**
  9. **[run_bcf_sensitivity_n100_dgp_c.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_bcf_sensitivity_n100_dgp_c.ipynb)**
  10. **[run_bcf_sensitivity_n250_dgp_a.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_bcf_sensitivity_n250_dgp_a.ipynb)**
  11. **[run_bcf_sensitivity_n250_dgp_b.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_bcf_sensitivity_n250_dgp_b.ipynb)**
  12. **[run_bcf_sensitivity_n250_dgp_c.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_bcf_sensitivity_n250_dgp_c.ipynb)**
  13. **[run_bcf_sensitivity_n1000_dgp_a.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_bcf_sensitivity_n1000_dgp_a.ipynb)**
  14. **[run_bcf_sensitivity_n1000_dgp_b.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_bcf_sensitivity_n1000_dgp_b.ipynb)**
  15. **[run_bcf_sensitivity_n1000_dgp_c.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_bcf_sensitivity_n1000_dgp_c.ipynb)**
- **ZIC-BCF-Smear**:
  16. **[run_zicbcf_smear_sensitivity_n100_dgp_a.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_zicbcf_smear_sensitivity_n100_dgp_a.ipynb)**
  17. **[run_zicbcf_smear_sensitivity_n100_dgp_b.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_zicbcf_smear_sensitivity_n100_dgp_b.ipynb)**
  18. **[run_zicbcf_smear_sensitivity_n100_dgp_c.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_zicbcf_smear_sensitivity_n100_dgp_c.ipynb)**
  19. **[run_zicbcf_smear_sensitivity_n250_dgp_a.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_zicbcf_smear_sensitivity_n250_dgp_a.ipynb)**
  20. **[run_zicbcf_smear_sensitivity_n250_dgp_b.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_zicbcf_smear_sensitivity_n250_dgp_b.ipynb)**
  21. **[run_zicbcf_smear_sensitivity_n250_dgp_c.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_zicbcf_smear_sensitivity_n250_dgp_c.ipynb)**
  22. **[run_zicbcf_smear_sensitivity_n1000_dgp_a.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_zicbcf_smear_sensitivity_n1000_dgp_a.ipynb)**
  23. **[run_zicbcf_smear_sensitivity_n1000_dgp_b.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_zicbcf_smear_sensitivity_n1000_dgp_b.ipynb)**
  24. **[run_zicbcf_smear_sensitivity_n1000_dgp_c.ipynb](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/colab_notebooks/run_zicbcf_smear_sensitivity_n1000_dgp_c.ipynb)**

### 2.3 Zero-Inflation (ZI) Sensitivity Analysis ($N=500$, non-standard $c_{\text{shift}}$ Levels 1, 2, 4, 5)
- **Standard BCF (Linear)**:
  25-36. `run_bcf_sensitivity_zi_lvl[1,2,4,5]_dgp_[a,b,c].ipynb` (12 notebooks)
- **ZIC-BCF-Smear**:
  37-48. `run_zicbcf_smear_sensitivity_zi_lvl[1,2,4,5]_dgp_[a,b,c].ipynb` (12 notebooks)

---

## 3. Local Simulation Framework

We updated the local R script environment to match the exact mathematical specifications and MCMC settings used in Colab, ensuring perfect comparability between local executions and Colab-generated CSVs:

* **[run_updated_simulation_study.R](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/run_updated_simulation_study.R)**: Master runner that evaluates Standard BCF and ZIC-BCF-Smear across the 100 seeds for the 3 main DGPs.
* **[run_updated_sensitivity.R](file:///home/hugo_souto/Stuff/Research/ZI-BCF/simulation_studies/run_updated_sensitivity.R)**: Local runner for N-sensitivity and ZI-sensitivity, saving checkpoints sequentially.

### How to Run Locally

You can execute the master simulation study locally using:
```bash
Rscript simulation_studies/run_updated_simulation_study.R
```
And execute the local sensitivity runner using:
```bash
Rscript simulation_studies/run_updated_sensitivity.R
```

All results will be saved as standard CSV files under the `simulation_studies/results/` directory. These CSVs share the exact same column names and metrics as the Colab-generated CSVs, allowing you to easily process them together using your local plotting and summary R scripts.
