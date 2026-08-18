# Independent-seed diagnostic for the Gamma-hurdle benchmark used in
# run_nonoverlap_calendar_day_analysis.R. That analysis now runs four chains
# with seeds 1 to 4 and reports formal convergence diagnostics; this script is
# retained as an out-of-sample check on a seed the main analysis never uses.

suppressPackageStartupMessages({
  library(zicbcf)
  library(readr)
})

script_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))])
if (length(script_file) != 1L) stop("Run this script with Rscript.")

data_dir <- dirname(normalizePath(script_file))
source(file.path(data_dir, "..", "oncology_common.R"))
out_dir <- file.path(data_dir, "analysis")
df <- read_csv(file.path(out_dir, "nonoverlap_calendar_day_analysis_data.csv"), show_col_types = FALSE)

y <- df$nonoverlap_severe_ae_days
z <- as.integer(df$treatment)
X <- onc_design_matrix(df)

set.seed(2026)
fit <- gamma_hurdle(y = y, z = z, x = X, nburn = 1000, nsim = 1000, nthin = 1)
hurdle_ate_draws <- rowMeans(fit$p1 - fit$p0)

summarise_draws <- function(draws) {
  interval <- quantile(draws, c(0.025, 0.975))
  c(Estimate = mean(draws), CI_low = interval[1L], CI_high = interval[2L], P_gt_0 = mean(draws > 0))
}

results <- rbind(
  c(Estimand = "Total non-overlapping severe-AE days", summarise_draws(fit$ate)),
  c(Estimand = "Pr(any non-overlapping severe-AE day)", summarise_draws(hurdle_ate_draws))
)
write_csv(as.data.frame(results), file.path(out_dir, "nonoverlap_gamma_seed_sensitivity.csv"))
print(results)
