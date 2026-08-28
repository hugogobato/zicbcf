## Exp B2/C/D: (i) fixed diagnostics incl. proper lag-1; (ii) multi-study
## replication with fresh seeds; (iii) extended-burn-in causal test.
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

fit_one <- function(seed, nburn = 2000, nsim = 4000) {
  set.seed(seed)
  t0 <- Sys.time()
  f <- zicbcf_smear(y = y, z = z, x_control = X, x_moderate = X,
                    pihat = pihat, nburn = nburn, nsim = nsim)
  list(sigma = as.numeric(f$sigma_c),
       mu_c_mean = rowMeans(f$mu_c),
       mins = round(as.numeric(Sys.time() - t0, units = "mins"), 2))
}

gw <- function(x, f1 = 0.1, f2 = 0.5) {
  n <- length(x)
  w1 <- x[seq_len(floor(f1 * n))]; w2 <- x[(floor((1 - f2) * n) + 1):n]
  se <- sqrt(var(w1)/length(w1) + var(w2)/length(w2))
  # iid SEs are appropriate here (rho_1 ~ .04)
  (mean(w1) - mean(w2)) / se
}
early_shift <- function(x) {  # standardized (first 10% - rest) / se of difference
  n <- length(x); k <- floor(0.1 * n)
  se <- sqrt(var(x[seq_len(k)]) / k + var(x[-seq_len(k)]) / (n - k))
  (mean(x[seq_len(k)]) - mean(x[-seq_len(k)])) / se
}
rho1 <- function(x) as.numeric(acf(x, lag.max = 1, plot = FALSE)$acf[2])

mode <- Sys.getenv("MODE", "replicate")

if (mode == "analyze") {
  chains <- lapply(seq_len(4), function(i)
    readRDS(sprintf("/tmp/opencode/geweke/onc_chain%d_sigma.rds", i)))
  tab <- NULL
  for (i in seq_along(chains)) {
    s <- chains[[i]]$sigma_c; m <- chains[[i]]$mu_c_mean
    tab <- rbind(tab, data.frame(
      quantity = "sigma_c", chain = i,
      z = gw(s), early = early_shift(s), rho_1 = rho1(s),
      level = mean(s)))
    tab <- rbind(tab, data.frame(
      quantity = "mean_mu_c", chain = i,
      z = gw(m), early = early_shift(m), rho_1 = rho1(m),
      level = mean(m)))
  }
  print(tab, digits = 4, row.names = FALSE)
} else {
  stid <- Sys.getenv("STUDY", "1"); nburn <- as.integer(Sys.getenv("NBURN", "2000"))
  seeds <- if (stid == "ext") 1:4 else if (stid %in% c("r1","r2","r3","r4","r5"))
    as.integer(seq(as.integer(gsub("\\D","",stid)) * 10 + 1,
                   length.out = 4))
  else if (stid == "base") 1:4 else {
    raw <- Sys.getenv("SEEDS", "41,42,43,44")
    if (grepl(",", raw)) as.integer(strsplit(raw, ",")[[1]]) else
      do.call(seq, as.list(as.integer(strsplit(raw, ":")[[1]])))
  }
  res <- NULL
  for (sd in seeds) {
    f <- fit_one(sd, nburn = nburn)
    res <- rbind(res, data.frame(
      study = stid, nburn = nburn, seed = sd,
      sigma_z = gw(f$sigma), sigma_early = early_shift(f$sigma),
      sig_level = mean(f$sigma), sig_sd = sd(f$sigma),
      muc_z = gw(f$mu_c_mean), muc_early = early_shift(f$mu_c_mean),
      mins = f$mins))
    cat("done", stid, "nburn", nburn, "seed", sd,
        sprintf("sigma_z=%.2f early=%+.2f | muc_z=%.2f early=%+.2f\n",
                tail(res$sigma_z,1), tail(res$sigma_early,1),
                tail(res$muc_z,1), tail(res$muc_early,1)))
    saveRDS(list(sigma = f$sigma, mu_c_mean = f$mu_c_mean),
            sprintf("/tmp/opencode/geweke/study_%s_nburn%d_seed%d.rds", stid, nburn, sd))
  }
  saveRDS(res, sprintf("/tmp/opencode/geweke/study_%s_nburn%d_summary.rds", stid, nburn))
  print(res, digits = 3, row.names = FALSE)
}
