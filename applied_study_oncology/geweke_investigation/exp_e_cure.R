## Cure-check for the mu_c-flagged chains + decay profile of flagged sigma chain
suppressMessages(library(coda))
setwd("/home/hugo_souto/Stuff/Research/zicbcf")
suppressMessages(library(zicbcf))
outdir <- "applied_study_oncology"
source(file.path(outdir, "oncology_common.R"), chdir = TRUE)
df <- read.csv(file.path(outdir, "zic_bcf_headneck_analysis_data.csv"),
               stringsAsFactors = FALSE)
y <- df$cumulative_severe_ae_duration
z <- as.integer(df$treatment)
X <- onc_design_matrix(df)
pihat <- rep(0.5, length(z))

gw <- function(x, f1 = 0.1, f2 = 0.5) {
  n <- length(x)
  w1 <- x[seq_len(floor(f1 * n))]; w2 <- x[(floor((1 - f2) * n) + 1):n]
  se <- sqrt(var(w1)/length(w1) + var(w2)/length(w2))
  (mean(w1) - mean(w2)) / se
}

for (sd_ in c(11, 12)) {
  set.seed(sd_)
  f <- zicbcf_smear(y = y, z = z, x_control = X, x_moderate = X,
                    pihat = pihat, nburn = 5000, nsim = 4000)
  m <- rowMeans(f$mu_c)
  saveRDS(list(sigma = as.numeric(f$sigma_c), mu_c_mean = m),
          sprintf("/tmp/opencode/geweke/cure_seed%d_nburn5000.rds", sd_))
  cat(sprintf("seed %d @nburn=5000: muc_z=%.2f (was -3.59 / -3.51 @2000)\n",
              sd_, gw(m)))
}

## Decay profile of the originally flagged chain (base seed 1, nburn=2000)
b <- readRDS("/tmp/opencode/geweke/study_base_nburn2000_seed1.rds")
s <- b$sigma
brk <- seq(0, length(s), length.out = 41)   # 2.5% blocks
wmean <- tapply(s, cut(seq_along(s), brk, include.lowest = TRUE), mean)
gs <- length(s) %/% 20
cat("\nsigma_c level (sd units above global mean) by 2.5%-block:\n")
print(round((as.numeric(wmean) - mean(s)) / (sd(s)/sqrt(gs)), 2))
