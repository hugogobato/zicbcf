# DESCRIPTIVE person-time analysis of the severe-adverse-event endpoints.
#
# Read the whole of this header before reading any number this script prints.
#
# WHAT THIS IS. The endpoint is rescaled to distinct severe-adverse-event days
# per 100 treatment-days, that is, the accrued burden divided by the patient's
# own time on treatment. The randomized analyses in run_zicbcf_oncology.R and
# run_nonoverlap_calendar_day_analysis.R are left untouched; this is reported
# alongside them.
#
# WHAT THIS IS NOT. It is not an adjusted causal effect and it must not be
# reported as one. Treatment duration is measured after randomization and is a
# consequence of assignment, so dividing by it conditions on a post-treatment
# variable and can induce selection bias. The direction of that bias is not
# fixed in advance and cannot be recovered from the supplied files. A patient
# who stops early because of toxicity contributes a short denominator alongside
# a large numerator, which inflates that patient's rate and makes their arm look
# worse; a patient who stops early for unrelated reasons pushes the same arm's
# rate down. Which of the two dominates depends on why patients discontinued,
# and the supplied extract carries no discontinuation-reason field.
#
# The two normalizations also disagree, which is itself the finding. The
# patient-level mean of per-patient rates weights every patient equally and so
# up-weights short-exposure patients; the aggregate rate pools days over pooled
# treatment-time. The defensible summary of this script is severe attenuation
# of the arm difference once exposure is accounted for, not reversal.
#
# A causal quantity that removes the duration channel from these data would
# need an explicit treatment-duration policy, stated as an assumption. An
# offset is not such an assumption.

suppressPackageStartupMessages({
  library(zicbcf)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

script_file <- sub("^--file=", "",
                   commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))])
if (length(script_file) != 1L) stop("Run this script with Rscript.")
data_dir <- dirname(normalizePath(script_file))
source(file.path(data_dir, "..", "oncology_common.R"))
out_dir <- file.path(data_dir, "analysis")
data_file <- file.path(out_dir, "nonoverlap_calendar_day_analysis_data.csv")
if (!file.exists(data_file)) stop("Run build_nonoverlap_calendar_day_endpoint.R first.")

df <- read_csv(data_file, show_col_types = FALSE)
if (!"treatment_duration" %in% names(df)) {
  stop("Rebuild the endpoint file: treatment_duration is required.")
}
if (any(df$treatment_duration <= 0)) stop("Treatment duration must be positive.")

z <- as.integer(df$treatment)
X <- onc_design_matrix(df)
labels <- onc_subgroup_labels(df)

N_CHAINS <- 4L
CHAIN_SEEDS <- seq_len(N_CHAINS)

# ---------------------------------------------------------------------------
# 1. Unadjusted description under both normalizations
# ---------------------------------------------------------------------------
descriptive <- df %>%
  mutate(rate_per_100 = 100 * nonoverlap_severe_ae_days / treatment_duration) %>%
  group_by(Arm = treatment_label) %>%
  summarise(
    n = n(),
    `Mean treatment days` = mean(treatment_duration),
    `Mean distinct severe-AE days` = mean(nonoverlap_severe_ae_days),
    `Patient-level mean rate per 100 treatment-days` = mean(rate_per_100),
    `Aggregate rate per 100 treatment-days` =
      100 * sum(nonoverlap_severe_ae_days) / sum(treatment_duration),
    .groups = "drop"
  )
write_csv(descriptive, file.path(out_dir, "exposure_rate_descriptive_summary.csv"))

# ---------------------------------------------------------------------------
# 2. The same two-part models, fit to the per-exposure rate
#
# The rate inherits the semicontinuous structure of the day-scale endpoint: the
# same patients contribute an exact zero, and the positive part remains heavily
# right-skewed. Fitting the same models to it isolates the effect of the change
# of scale from the effect of the change of estimator.
# ---------------------------------------------------------------------------
y_rate <- 100 * df$nonoverlap_severe_ae_days / df$treatment_duration

cat(sprintf("Rate outcome: %d zeros of %d (%.1f%%); positive median %.2f, max %.2f\n",
            sum(y_rate == 0), length(y_rate), 100 * mean(y_rate == 0),
            median(y_rate[y_rate > 0]), max(y_rate)))

pihat <- rep(0.5, length(z))

zic_chains <- lapply(CHAIN_SEEDS, function(seed) {
  set.seed(seed)
  cat("  ZIC-BCF-Smear chain", seed, "\n")
  zicbcf_smear(y = y_rate, z = z, x_control = X, x_moderate = X, pihat = pihat,
               nburn = 2000, nsim = 4000)
})
zic_convergence <- zicbcf_convergence(zic_chains, n_cate_units = 5L)
zic_ate <- unlist(lapply(zic_chains, function(f) f$ate), use.names = FALSE)
zic_cate <- do.call(rbind, lapply(zic_chains, function(f) f$cate))

gamma_chains <- lapply(CHAIN_SEEDS, function(seed) {
  set.seed(seed)
  cat("  Gamma hurdle chain", seed, "\n")
  gamma_hurdle(y = y_rate, z = z, x = X, nburn = 1000, nsim = 1000, nthin = 1)
})
gamma_ate <- unlist(lapply(gamma_chains, function(f) f$ate), use.names = FALSE)
gamma_cate <- do.call(rbind, lapply(gamma_chains, function(f) f$cate))

results <- bind_rows(
  data.frame(Model = "ZIC-BCF-Smear",
             Estimand = "DESCRIPTIVE rate: severe-AE days per 100 treatment-days",
             t(onc_posterior_summary(zic_ate))),
  data.frame(Model = "Gamma hurdle",
             Estimand = "DESCRIPTIVE rate: severe-AE days per 100 treatment-days",
             t(onc_posterior_summary(gamma_ate)))
)
write_csv(results, file.path(out_dir, "exposure_rate_model_results.csv"))
write_csv(zic_convergence, file.path(out_dir, "exposure_rate_convergence_diagnostics.csv"))

subgroups <- bind_rows(
  onc_subgroup_table(zic_cate, labels, "ZIC-BCF-Smear",
                     "DESCRIPTIVE rate per 100 treatment-days"),
  onc_subgroup_table(gamma_cate, labels, "Gamma hurdle",
                     "DESCRIPTIVE rate per 100 treatment-days")
)
write_csv(subgroups, file.path(out_dir, "exposure_rate_subgroup_results.csv"))

saveRDS(list(
  outcome = "Distinct severe-AE days per 100 treatment-days (DESCRIPTIVE)",
  warning = paste("Treatment duration is post-randomization. These are not",
                  "adjusted causal effects and must not be reported as such."),
  zic = list(ate = zic_ate, cate = zic_cate),
  gamma = list(ate = gamma_ate, cate = gamma_cate)
), file.path(out_dir, "exposure_rate_posterior_draws.rds"))

rate_plot_data <- df %>%
  mutate(rate_per_100 = 100 * nonoverlap_severe_ae_days / treatment_duration)
png(file.path(out_dir, "exposure_rate_distribution.png"), width = 800, height = 550)
print(
  ggplot(rate_plot_data, aes(x = treatment_label, y = rate_per_100,
                             fill = treatment_label)) +
    geom_boxplot(outlier.alpha = 0.3) +
    scale_y_continuous(trans = "sqrt") +
    labs(x = NULL, y = "Distinct severe-AE days per 100 treatment-days (square-root scale)",
         title = "Descriptive per-exposure burden by randomized arm",
         subtitle = "Descriptive only: treatment duration is post-randomization") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
)
dev.off()

cat("\n=== Descriptive summary by arm ===\n")
print(as.data.frame(descriptive), row.names = FALSE, digits = 4)
cat("\n=== Two-part models fit to the per-exposure rate (DESCRIPTIVE) ===\n")
print(as.data.frame(results), row.names = FALSE, digits = 4)
cat(sprintf("\nZIC-BCF-Smear: max Rhat = %.4f, min ESS = %.0f\n",
            max(zic_convergence$rhat, na.rm = TRUE),
            min(zic_convergence$ess, na.rm = TRUE)))
cat("\nReminder: these are descriptive rates per treatment-day, not adjusted causal effects.\n")
