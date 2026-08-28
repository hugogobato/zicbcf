## Exp B analysis: dissect the four real sigma_c chains
suppressMessages(library(coda))
set.seed(20260826)

geweke_z <- function(x, f1 = 0.1, f2 = 0.5) {
  as.numeric(coda::geweke.diag(coda::mcmc(as.numeric(x)), frac1 = f1, frac2 = f2)$z)
}

chains <- lapply(seq_len(4), function(i)
  readRDS(sprintf("/tmp/opencode/geweke/onc_chain%d_sigma.rds", i)))

pdf("/tmp/opencode/geweke/sigma_c_dissection.pdf", width = 10, height = 8)
par(mfrow = c(3, 4), mar = c(3, 3, 2, 1))
all_z <- list(); z_table <- NULL
for (i in seq_len(4)) {
  s <- chains[[i]]$sigma_c; n <- length(s)

  # trace plot
  plot(s, type = "l", main = sprintf("chain %d sigma_c (n=%d)", i, n),
       xlab = "draw", ylab = "sigma_c")

  # ACF: negative autocorrelation signature?
  a <- acf(s, lag.max = 40, plot = TRUE)$acf

  # rolling means over successive deciles: drift or noise?
  dec <- cut(seq_along(s), breaks = seq(0, n, length.out = 21), include.lowest = TRUE)
  mw <- tapply(s, dec, mean)
  if (i == 1) {
    cat("\n--- Within-chain drift check: mean of each successive 5% block ---\n")
    print(round(mw, 5)); cat("sd of block means:", sd(mw), "\n")
  }

  # Geweke across window choices
  zs <- c(
    std        = geweke_z(s, .1, .5),
    balanced   = geweke_z(s, .25, .25),
    lasthalf   = NA,
    f02_f08    = geweke_z(s, .02, .8),
    iidSE_01_05 = (mean(s[seq_len(floor(.1 * n))]) -
                   mean(s[(floor(.5 * n) + 1):n])) /
                  (sd(s) * sqrt(1 / floor(.1 * n) + 1 / (n - floor(.5 * n))))
  )
  zs["lasthalf"] <- (mean(s[seq_len(floor(.1*n))]) - mean(s)) / (sd(s) * sqrt((1/floor(.1*n))))
  all_z[[i]] <- zs

  # spectral density at zero estimate vs naive variance/n: the ratio is what
  # gweke's denominator uses per window. Report var est for first-10% window.
  w1 <- s[seq_len(floor(.1 * n))]
  w2 <- s[(floor(.5 * n) + 1):n]
  s0_w1 <- coda::spectrum0.ar(coda::mcmc(w1))$spec
  s0_w2 <- coda::spectrum0.ar(coda::mcmc(w2))$spec
  se_spec <- sqrt(s0_w1 / length(w1) + s0_w2 / length(w2))
  se_naive <- sqrt(var(w1) / length(w1) + var(w2) / length(w2))

  # effectiveSize per chain
  ess <- coda::effectiveSize(coda::mcmc(s))
  ac1 <- mean(acf(s, lag.max = 1, plot = FALSE)$acf)

  z_table <- rbind(z_table, data.frame(
    chain = i, z_std = zs["std"], z_balanced = zs["balanced"],
    z_f02f08 = zs["f02_f08"], z_iidSE = zs["iidSE_01_05"],
    ess_ratio = ess / n, lag1_acf = ac1,
    spec_over_var_w1 = s0_w1 / var(w1), spec_over_var_w2 = s0_w2 / var(w2),
    mean_first10 = mean(w1), sd_first10 = sd(w1), sd_last50 = sd(w2)))
}
dev.off()

cat("\n=== per-chain summary ===\n"); print(z_table, row.names = FALSE)
cat("\nz under null max|z| of 16 chain-level stats:",
    max(abs(unlist(all_z))), "\n")

## Also: does ANY quantity pattern replicate? compare our chains' full monitored set
library(zicbcf)
