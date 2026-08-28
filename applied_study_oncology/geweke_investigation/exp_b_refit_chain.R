## Exp B: reproduce the oncology fit and dissect the sigma_c chains.
## Same data, same config as run_zicbcf_oncology.R (4 chains, nburn=2000,
## nsim=4000), but we KEEP the per-chain sigma_c draws for analysis.
setwd("/home/hugo_souto/Stuff/Research/zicbcf")
suppressMessages(library(zicbcf))
outdir <- "applied_study_oncology"
source(file.path(outdir, "oncology_common.R"))

df <- read.csv(file.path(outdir, "zic_bcf_headneck_analysis_data.csv"),
               stringsAsFactors = FALSE)
y <- df$cumulative_severe_ae_duration
z <- as.integer(df$treatment)
X <- onc_design_matrix(df)
pihat <- rep(0.5, length(z))

chain_id <- as.integer(Sys.getenv("CHAIN_ID", "1"))
seed <- chain_id  # CHAIN_SEEDS <- seq_len(N_CHAINS)
set.seed(seed)

t0 <- Sys.time()
fit <- zicbcf_smear(y = y, z = z, x_control = X, x_moderate = X,
                    pihat = pihat, nburn = 2000, nsim = 4000)
cat("elapsed:", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")

saveRDS(list(sigma_c = as.numeric(fit$sigma_c),
             smearing_factors = fit$smearing_factors,
             ate = fit$ate,
             mu_c_mean = rowMeans(fit$mu_c),
             tau_c_mean = rowMeans(fit$tau_c),
             seed = seed),
        sprintf("/tmp/opencode/geweke/onc_chain%d_sigma.rds", chain_id))
cat("saved chain", chain_id, "\n")
