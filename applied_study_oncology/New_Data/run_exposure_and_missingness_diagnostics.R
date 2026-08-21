# Descriptive diagnostics for the SPECTRUM severe-adverse-event endpoints.
#
# Nothing in this script is a causal estimate, and nothing in it is intended to
# correct one. Its purpose is to make three properties of the data visible that
# the model fits cannot express.
#
# 1. TREATMENT EXPOSURE. The endpoints count days on which a patient was in a
#    severe-adverse-event state. That quantity accrues with time on treatment.
#    Panitumumab continues as maintenance after chemotherapy stops, so time on
#    treatment is structurally longer in the experimental arm. Treatment
#    duration is a post-randomization consequence of assignment: it lies on the
#    causal pathway from assignment to accrued burden, so it can be neither
#    ignored nor naively conditioned on. What this script reports is how much of
#    the observed arm difference is exposure-dependent.
#
#    The common-horizon table truncates each patient's observed accrual at a
#    common study day. It is a diagnostic, not a fixed-horizon estimand.
#    Truncating at day 84 removes events recorded after day 84, but it does not
#    create observation for a patient whose follow-up ended before then, and the
#    supplied files carry no patient-specific follow-up end date. Read the table
#    as evidence that the effect is exposure-dependent, never as a corrected
#    estimate.
#
# 2. AN UNRESOLVED INTERVAL BOUNDARY. The distinct-calendar-day endpoint keeps
#    only qualifying records with a complete imputed interval and assigns zero
#    to patients with none. This script reports exactly how many records and
#    patients that affects, and how the affected patients split across arms.
#
# 3. PRIOR TREATMENT FOR HEAD-AND-NECK CANCER. PRHNTRTC is measured, prevalent
#    and imbalanced, and prior therapy to the head and neck is a principal
#    determinant of the mucosal and cutaneous toxicity being modeled. Its
#    distribution is reported here. The v2 domain diagnostics separately
#    report the now-available radiotherapy-specific PRRADIO variable used in
#    oncology_common.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

script_file <- sub("^--file=", "",
                   commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))])
if (length(script_file) != 1L) stop("Run this script with Rscript.")
data_dir <- dirname(normalizePath(script_file))
out_dir <- file.path(data_dir, "analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ae_path <- file.path(data_dir, "ae.sas7bdat")
patient_path <- file.path(data_dir, "..", "zic_bcf_headneck_analysis_data.csv")
raw_patient_path <- file.path(data_dir, "..", "headneck_265_severe_ae_patient_level.csv")
if (!file.exists(ae_path)) stop("Missing AE-level file: ", ae_path)
if (!file.exists(patient_path)) stop("Run prepare_analysis_data.R first.")

ae <- read_sas(ae_path) %>%
  mutate(across(where(haven::is.labelled), haven::zap_labels),
         SUBJID = as.character(SUBJID))
raw_patients <- read_csv(raw_patient_path, show_col_types = FALSE) %>%
  mutate(SUBJID = as.character(SUBJID),
         arm = if_else(TRTCD == 2, "Panitumumab + chemotherapy", "Chemotherapy"))
patients <- read_csv(patient_path, show_col_types = FALSE)

ARMS <- c("Chemotherapy", "Panitumumab + chemotherapy")

# ---------------------------------------------------------------------------
# Union of inclusive day intervals, optionally truncated at a common horizon.
# ---------------------------------------------------------------------------
union_days <- function(start_day, end_day, horizon = Inf) {
  start_day <- pmax(start_day, 1)
  end_day <- pmin(end_day, horizon)
  keep <- start_day <= end_day
  start_day <- start_day[keep]
  end_day <- end_day[keep]
  if (!length(start_day)) return(0)
  order_index <- order(start_day, end_day)
  start_day <- start_day[order_index]
  end_day <- end_day[order_index]
  total <- 0
  current_start <- start_day[1L]
  current_end <- end_day[1L]
  if (length(start_day) > 1L) {
    for (index in seq.int(2L, length(start_day))) {
      if (start_day[index] <= current_end + 1L) {
        current_end <- max(current_end, end_day[index])
      } else {
        total <- total + current_end - current_start + 1L
        current_start <- start_day[index]
        current_end <- end_day[index]
      }
    }
  }
  total + current_end - current_start + 1L
}

qualifying <- ae %>%
  mutate(severe = AESEVCD >= 3,
         on_study = is.na(AEPRIOR) | AEPRIOR != "Y") %>%
  filter(severe, on_study)

# ---------------------------------------------------------------------------
# 1. Cohort reconciliation
# ---------------------------------------------------------------------------
cohort <- tibble(
  Quantity = c(
    "Randomized in the parent trial (Vermorken et al. 2013)",
    "Patients in the analytic cohort",
    "Analytic cohort, chemotherapy",
    "Analytic cohort, panitumumab plus chemotherapy",
    "Distinct subjects in the adverse-event file",
    "Adverse-event file subjects, chemotherapy",
    "Adverse-event file subjects, panitumumab plus chemotherapy",
    "Analytic cohort patients with no adverse-event record"
  ),
  Value = c(
    657,
    nrow(patients),
    sum(patients$treatment == 0L),
    sum(patients$treatment == 1L),
    n_distinct(ae$SUBJID),
    sum(!duplicated(ae$SUBJID) & ae$TRTCD == 50),
    sum(!duplicated(ae$SUBJID) & ae$TRTCD == 2),
    nrow(patients) - n_distinct(ae$SUBJID)
  )
)
write_csv(cohort, file.path(out_dir, "cohort_reconciliation.csv"))

# ---------------------------------------------------------------------------
# 2. Treatment exposure
# ---------------------------------------------------------------------------
exposure <- patients %>%
  group_by(Arm = treatment_label) %>%
  summarise(
    n = n(),
    Mean = mean(treatment_duration),
    SD = sd(treatment_duration),
    Median = median(treatment_duration),
    Q1 = quantile(treatment_duration, 0.25),
    Q3 = quantile(treatment_duration, 0.75),
    Minimum = min(treatment_duration),
    Maximum = max(treatment_duration),
    .groups = "drop"
  )
exposure_ratio <- with(exposure, Mean[Arm == ARMS[2]] / Mean[Arm == ARMS[1]])
write_csv(exposure, file.path(out_dir, "treatment_exposure_by_arm.csv"))

# ---------------------------------------------------------------------------
# 3. Common-horizon truncation table, across all four start/end date variants
# ---------------------------------------------------------------------------
variants <- list(
  `Imputed start, imputed end (primary)` = c("AESTDYI", "AEENDYI"),
  `Raw start, raw end` = c("AESTDY", "AEENDY"),
  `Raw start, imputed end` = c("AESTDY", "AEENDYI"),
  `Imputed start, raw end` = c("AESTDYI", "AEENDY")
)
horizons <- c(Inf, 120, 90, 84)

truncation <- bind_rows(lapply(names(variants), function(variant_name) {
  columns <- variants[[variant_name]]
  complete <- qualifying %>%
    filter(is.finite(.data[[columns[1L]]]), is.finite(.data[[columns[2L]]]))
  bind_rows(lapply(horizons, function(horizon) {
    per_patient <- complete %>%
      group_by(SUBJID) %>%
      summarise(days = union_days(.data[[columns[1L]]], .data[[columns[2L]]], horizon),
                .groups = "drop")
    merged <- raw_patients %>%
      select(SUBJID, arm) %>%
      left_join(per_patient, by = "SUBJID") %>%
      mutate(days = coalesce(days, 0))
    means <- merged %>% group_by(arm) %>%
      summarise(mean_days = mean(days), .groups = "drop")
    tibble(
      Variant = variant_name,
      Horizon = if (is.infinite(horizon)) "Uncapped" else paste(horizon, "days"),
      Chemotherapy = means$mean_days[means$arm == ARMS[1]],
      Panitumumab = means$mean_days[means$arm == ARMS[2]],
      Difference = means$mean_days[means$arm == ARMS[2]] -
        means$mean_days[means$arm == ARMS[1]]
    )
  }))
}))
write_csv(truncation, file.path(out_dir, "common_horizon_truncation.csv"))

# ---------------------------------------------------------------------------
# 4. Person-time rates
#
# Two normalizations are reported because they do not agree, and the
# disagreement is the point. The patient-level mean of per-patient rates gives
# every patient equal weight and therefore up-weights short-exposure patients;
# the aggregate rate pools days over pooled treatment-time. Neither is an
# adjusted causal effect. Dividing by treatment duration conditions on a
# post-randomization quantity, which can induce selection bias whose direction
# depends on why patients discontinued, and that is not recoverable from the
# supplied files.
# ---------------------------------------------------------------------------
primary_days <- qualifying %>%
  filter(is.finite(AESTDYI), is.finite(AEENDYI)) %>%
  group_by(SUBJID) %>%
  summarise(distinct_days = union_days(AESTDYI, AEENDYI), .groups = "drop")

per_patient <- raw_patients %>%
  select(SUBJID, arm, TRTDUR, severe_ae_total_duration) %>%
  left_join(primary_days, by = "SUBJID") %>%
  mutate(distinct_days = coalesce(distinct_days, 0),
         rate_per_100 = 100 * distinct_days / TRTDUR)

person_time <- per_patient %>%
  group_by(Arm = arm) %>%
  summarise(
    n = n(),
    `Mean treatment days` = mean(TRTDUR),
    `Total treatment days` = sum(TRTDUR),
    `Total distinct severe-AE days` = sum(distinct_days),
    `Patient-level mean rate per 100 treatment-days` = mean(rate_per_100),
    `Aggregate rate per 100 treatment-days` = 100 * sum(distinct_days) / sum(TRTDUR),
    .groups = "drop"
  )
write_csv(person_time, file.path(out_dir, "person_time_rates.csv"))

# ---------------------------------------------------------------------------
# 5. Unresolved interval boundaries
# ---------------------------------------------------------------------------
qualifying <- qualifying %>%
  mutate(complete_interval = is.finite(AESTDYI) & is.finite(AEENDYI))
incomplete <- qualifying %>% filter(!complete_interval)

by_patient <- qualifying %>%
  group_by(SUBJID) %>%
  summarise(qualifying_records = n(),
            complete_records = sum(complete_interval), .groups = "drop") %>%
  left_join(raw_patients %>% select(SUBJID, arm), by = "SUBJID")

forced_zero <- by_patient %>% filter(complete_records == 0L)
affected <- by_patient %>% filter(complete_records < qualifying_records)

arm_count <- function(data, arm_name) sum(data$arm == arm_name)

boundary <- tibble(
  Quantity = c(
    "Qualifying on-study grade 3+ records",
    "Records with a complete imputed interval",
    "Records with an unresolved boundary",
    "Unresolved records missing the end day only",
    "Unresolved records missing the start day only",
    "Patients with at least one unresolved record",
    "  of whom chemotherapy",
    "  of whom panitumumab plus chemotherapy",
    "Patients with no complete qualifying interval, assigned zero",
    "  of whom chemotherapy",
    "  of whom panitumumab plus chemotherapy"
  ),
  Value = c(
    nrow(qualifying),
    sum(qualifying$complete_interval),
    nrow(incomplete),
    sum(is.finite(incomplete$AESTDYI) & !is.finite(incomplete$AEENDYI)),
    sum(!is.finite(incomplete$AESTDYI) & is.finite(incomplete$AEENDYI)),
    nrow(affected),
    arm_count(affected, ARMS[1]),
    arm_count(affected, ARMS[2]),
    nrow(forced_zero),
    arm_count(forced_zero, ARMS[1]),
    arm_count(forced_zero, ARMS[2])
  )
)
write_csv(boundary, file.path(out_dir, "unresolved_interval_boundaries.csv"))

# The file carries no ongoing-at-cutoff flag, so a missing end day cannot be
# shown to mean the event was still running when data collection stopped. It is
# not evidence about observation windows. What is measurable is the arm split.

# ---------------------------------------------------------------------------
# 6. Prior treatment for head-and-neck cancer
# ---------------------------------------------------------------------------
prior_prevalence <- patients %>%
  group_by(Arm = treatment_label) %>%
  summarise(
    n = n(),
    `Prior treatment (n)` = sum(prior_hn_treatment == "Y"),
    `Prior treatment (%)` = 100 * mean(prior_hn_treatment == "Y"),
    `IVRS prior treatment (n)` = sum(prior_hn_treatment_ivrs == "Y"),
    `IVRS prior treatment (%)` = 100 * mean(prior_hn_treatment_ivrs == "Y"),
    .groups = "drop"
  )
write_csv(prior_prevalence, file.path(out_dir, "prior_treatment_prevalence.csv"))

prior_concordance <- patients %>%
  count(prior_hn_treatment, prior_hn_treatment_ivrs) %>%
  rename(`PRHNTRTC` = prior_hn_treatment, `PRHNTRTI (IVRS)` = prior_hn_treatment_ivrs)
write_csv(prior_concordance, file.path(out_dir, "prior_treatment_concordance.csv"))

prior_stratified <- per_patient %>%
  left_join(raw_patients %>% select(SUBJID, PRHNTRTC), by = "SUBJID") %>%
  group_by(`Prior SCCHN treatment` = PRHNTRTC, Arm = arm) %>%
  summarise(
    n = n(),
    `Mean distinct severe-AE days` = mean(distinct_days),
    `Mean summed record days` = mean(severe_ae_total_duration),
    `Mean treatment days` = mean(TRTDUR),
    .groups = "drop"
  )
write_csv(prior_stratified, file.path(out_dir, "prior_treatment_stratified_burden.csv"))

# ---------------------------------------------------------------------------
# 7. Tumor category is a deterministic recode of diagnosis site
# ---------------------------------------------------------------------------
recode_map <- patients %>% count(diagtype, tumcat)
write_csv(recode_map, file.path(out_dir, "tumor_category_recode_map.csv"))

# ---------------------------------------------------------------------------
# 8. Figure: exposure distribution and the attenuation of the arm difference
# ---------------------------------------------------------------------------
primary_truncation <- truncation %>%
  filter(Variant == "Imputed start, imputed end (primary)") %>%
  mutate(Horizon = factor(Horizon, levels = c("84 days", "90 days", "120 days", "Uncapped")))

exposure_plot <- ggplot(patients, aes(x = treatment_label, y = treatment_duration,
                                      fill = treatment_label)) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_y_continuous(trans = "log10") +
  labs(x = NULL, y = "Treatment duration (days, log scale)",
       title = "Time on treatment by randomized arm") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

truncation_plot <- ggplot(primary_truncation,
                          aes(x = Horizon, y = Difference, group = 1)) +
  geom_line(colour = "firebrick", linewidth = 1) +
  geom_point(colour = "firebrick", size = 3) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey40") +
  labs(x = "Common truncation horizon", y = "Arm difference (mean distinct severe-AE days)",
       title = "Unadjusted arm difference under common-horizon truncation") +
  theme_minimal(base_size = 12)

png(file.path(out_dir, "treatment_exposure_distribution.png"), width = 750, height = 550)
print(exposure_plot)
dev.off()
png(file.path(out_dir, "common_horizon_attenuation.png"), width = 750, height = 550)
print(truncation_plot)
dev.off()

cat("=== Cohort reconciliation ===\n"); print(as.data.frame(cohort), row.names = FALSE)
cat("\n=== Treatment exposure by arm ===\n"); print(as.data.frame(exposure), row.names = FALSE)
cat(sprintf("Ratio of mean treatment duration (panitumumab / chemotherapy): %.3f\n",
            exposure_ratio))
cat("\n=== Common-horizon truncation ===\n")
print(as.data.frame(truncation), row.names = FALSE, digits = 4)
cat("\n=== Person-time rates ===\n"); print(as.data.frame(person_time), row.names = FALSE, digits = 4)
cat("\n=== Unresolved interval boundaries ===\n"); print(as.data.frame(boundary), row.names = FALSE)
cat("\n=== Prior SCCHN treatment ===\n"); print(as.data.frame(prior_prevalence), row.names = FALSE, digits = 4)
cat("\n=== Prior-treatment concordance (CRF composite versus IVRS) ===\n")
print(as.data.frame(prior_concordance), row.names = FALSE)
cat("\n=== Burden by prior treatment and arm ===\n")
print(as.data.frame(prior_stratified), row.names = FALSE, digits = 4)
cat("\n=== Tumor category by diagnosis site ===\n"); print(as.data.frame(recode_map), row.names = FALSE)
cat("\nDiagnostics completed.\n")
