## Exp A: null calibration of coda::geweke.diag(frac1=.1, frac2=.5)
## Question: under a perfectly stationary chain, what does the distribution of
## the Geweke z look like, and how anomalous is max|z| over 8 quantities x 4 chains?
library(coda)
set.seed(20260826)

geweke_z <- function(x) {
  suppressWarnings(as.numeric(coda::geweke.diag(coda::mcmc(as.numeric(x)),
                                                frac1 = 0.1, frac2 = 0.5)$z))
}

run_sim <- function(gen, n_draws, reps = 4000) {
  zs <- vapply(seq_len(reps), function(i) geweke_z(gen(n_draws)), numeric(1))
  zs[is.finite(zs)]
}

## A1: iid normal chains (the textbook H0)
n_reps <- 4000
res <- list()
for (nn in c(2000, 4000)) {
  z_iid <- run_sim(function(n) rnorm(n), nn, n_reps)
  res[[paste0("iid", nn)]] <- z_iid
  cat(sprintf("iid N(0,1), n=%d: mean=%.3f sd=%.3f  P(|z|>3.38)=%.4f  P(|z|>=3.67)=%.4f\n",
              nn, mean(z_iid), sd(z_iid),
              mean(abs(z_iid) > 3.38), mean(abs(z_iid) >= 3.67)))
}

## A2: antithetic AR(1) chains with negative rho (ESS > n case seen in the data:
## smearing factor ESS/n = 1.244 implies rho_1 ~ -0.1 to -0.2; sigma_c ESS/n ~ .9)
for (rho in c(-0.05, -0.1, -0.2, -0.35, -0.5)) {
  gen <- function(n) {
    e <- rnorm(n + 200); x <- as.numeric(stats::filter(e, rho, method = "recursive")); x[-(1:200)]
  }
  for (nn in c(2000, 4000)) {
    z_ar <- run_sim(gen, nn, n_reps)
    # theoretical ESS ratio for AR(1): (1+rho)/(1-rho)
    ess_ratio <- (1 + rho) / (1 - rho)
    cat(sprintf("AR(1) rho=%+.2f (ess/n~%.2f), n=%d: mean=%.3f sd=%.3f  P(|z|>3.38)=%.4f\n",
                rho, ess_ratio, nn, mean(z_ar), sd(z_ar), mean(abs(z_ar) > 3.38)))
  }
}

## A3: family-wise view -- max |z| over 72 *independent* near-iid deviates,
## i.e. exactly what the manuscripts' corrected critical value assumes.
zmax_null <- replicate(200000, max(abs(rnorm(72))))
cat(sprintf("\nNull max-of-72: P(max|z| >= 3.67) = %.4f (theory ~ %.4f)\n",
            mean(zmax_null >= 3.67), 1 - pnorm(3.67)^72))
cat(sprintf("Null max-of-72: P(max|z| >= 3.38) = %.4f\n", mean(zmax_null >= 3.38)))
cat(sprintf("median of max-of-72 = %.2f\n", median(zmax_null)))
