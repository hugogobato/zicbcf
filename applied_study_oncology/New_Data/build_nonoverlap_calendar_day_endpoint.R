# Construct the alternative severe-adverse-event endpoint from the AE-level
# SPECTRUM extract.  The estimand is the number of distinct study days on
# which a patient had at least one on-study grade 3+ AE.  Overlapping and
# adjacent intervals are counted once.  Dates are inclusive because the data
# dictionary defines AEDURI as (AEENDTI - AESTDTI) + 1.

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(readr)
  library(ggplot2)
})

script_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))])
if (length(script_file) != 1L) stop("Run this script with Rscript.")

data_dir <- dirname(normalizePath(script_file))
out_dir <- file.path(data_dir, "analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ae_path <- file.path(data_dir, "ae.sas7bdat")
patient_path <- file.path(data_dir, "..", "headneck_265_severe_ae_patient_level.csv")

if (!file.exists(ae_path)) stop("Missing AE-level file: ", ae_path)
if (!file.exists(patient_path)) stop("Missing patient-level file: ", patient_path)

# Return the number of integer study days covered by the union of inclusive
# intervals.  Adjacency is merged because [1, 3] and [4, 6] represent six
# uninterrupted calendar days with at least one qualifying AE.
union_days <- function(start_day, end_day) {
  stopifnot(length(start_day) == length(end_day), length(start_day) > 0L)
  order_index <- order(start_day, end_day)
  start_day <- start_day[order_index]
  end_day <- end_day[order_index]

  if (length(start_day) == 1L) return(end_day[1L] - start_day[1L] + 1L)

  total <- 0
  current_start <- start_day[1L]
  current_end <- end_day[1L]

  for (index in seq.int(2L, length(start_day))) {
    if (start_day[index] <= current_end + 1L) {
      current_end <- max(current_end, end_day[index])
    } else {
      total <- total + current_end - current_start + 1L
      current_start <- start_day[index]
      current_end <- end_day[index]
    }
  }

  total + current_end - current_start + 1L
}

ae <- read_sas(ae_path) %>%
  mutate(
    SUBJID = as.character(SUBJID),
    severe_ae = AESEVCD >= 3,
    on_study = is.na(AEPRIOR) | AEPRIOR != "Y",
    recorded_duration = coalesce(as.numeric(AEDUR), as.numeric(AEDURI)),
    complete_imputed_interval = is.finite(AESTDYI) & is.finite(AEENDYI)
  )

patients <- read_csv(patient_path, show_col_types = FALSE) %>%
  mutate(SUBJID = as.character(SUBJID))

if (nrow(patients) != n_distinct(patients$SUBJID)) {
  stop("The patient-level file must contain exactly one row per SUBJID.")
}

required_patient_variables <- c(
  "SUBJID", "AGE", "SEX", "TRTCD", "TRT", "B_ECOGCT", "DIAGTYPE",
  "TUMCAT", "HPV", "DSTATUS", "PRHNTRTC", "TRTDUR", "severe_ae_total_duration"
)
missing_patient_variables <- setdiff(required_patient_variables, names(patients))
if (length(missing_patient_variables) > 0L) {
  stop("Patient-level file is missing: ", paste(missing_patient_variables, collapse = ", "))
}

qualifying_events <- ae %>% filter(severe_ae, on_study)
complete_events <- qualifying_events %>% filter(complete_imputed_interval)

# Per-patient completeness of the qualifying records. Patients whose qualifying
# records all lack a resolvable interval are assigned zero distinct severe-AE
# days below; that assignment is a data limitation, not an observation, so the
# count is reported rather than left for a reader to infer.
patient_completeness <- qualifying_events %>%
  group_by(SUBJID) %>%
  summarise(qualifying_records = n(),
            complete_records = sum(complete_imputed_interval), .groups = "drop")

# First verify that the earlier record-duration endpoint is reproducible from
# this AE extract.  This also establishes a valid linkage to the baseline and
# randomization data in the patient-level file.
record_duration_by_patient <- qualifying_events %>%
  group_by(SUBJID) %>%
  summarise(
    recomputed_record_duration = sum(recorded_duration, na.rm = TRUE),
    .groups = "drop"
  )

record_duration_check <- patients %>%
  select(SUBJID, severe_ae_total_duration) %>%
  left_join(record_duration_by_patient, by = "SUBJID") %>%
  mutate(
    recomputed_record_duration = coalesce(recomputed_record_duration, 0),
    difference = severe_ae_total_duration - recomputed_record_duration
  )

if (any(record_duration_check$difference != 0)) {
  stop("The supplied AE file does not reproduce the prior record-duration endpoint.")
}

# For events with a documented imputed start and end study day, construct the
# patient-level union of the inclusive intervals.  An unresolved end date is
# not extended beyond the observed data because no patient-specific follow-up
# end date is available in the supplied extract.
calendar_days_by_patient <- complete_events %>%
  group_by(SUBJID) %>%
  summarise(
    complete_record_days = sum(AEENDYI - AESTDYI + 1),
    nonoverlap_severe_ae_days = union_days(AESTDYI, AEENDYI),
    overlap_days_removed = complete_record_days - nonoverlap_severe_ae_days,
    complete_severe_ae_records = n(),
    .groups = "drop"
  )

analysis_data <- patients %>%
  transmute(
    SUBJID,
    treatment = case_when(
      TRTCD == 2 ~ 1L,
      TRTCD == 50 ~ 0L,
      TRUE ~ NA_integer_
    ),
    treatment_label = case_when(
      TRTCD == 2 ~ "Panitumumab + chemotherapy",
      TRTCD == 50 ~ "Chemotherapy",
      TRUE ~ as.character(TRT)
    ),
    cumulative_severe_ae_duration = severe_ae_total_duration,
    age = AGE,
    sex = SEX,
    b_ecogct = B_ECOGCT,
    diagtype = DIAGTYPE,
    tumcat = TUMCAT,
    # The prior analytic file encodes absent baseline HPV status as "Unknown".
    # Retain that prespecified category so that the adjusted analysis has the
    # same 520-patient covariate set as the original record-duration analysis.
    hpv = coalesce(HPV, "Unknown"),
    dstatus = DSTATUS,
    # Prior treatment for squamous-cell carcinoma of the head and neck. Prior
    # therapy to the head and neck is a principal determinant of the mucosal
    # and cutaneous toxicity this endpoint measures, so it enters the
    # adjustment set (see ../oncology_common.R).
    prior_hn_treatment = PRHNTRTC,
    # Time on treatment. Post-randomization, therefore never adjusted for; it
    # is carried only so that the descriptive exposure analyses can run from
    # this file.
    treatment_duration = TRTDUR
  ) %>%
  left_join(calendar_days_by_patient, by = "SUBJID") %>%
  mutate(
    across(
      c(complete_record_days, nonoverlap_severe_ae_days,
        overlap_days_removed, complete_severe_ae_records),
      ~ coalesce(.x, 0)
    ),
    any_nonoverlap_severe_ae_day = as.integer(nonoverlap_severe_ae_days > 0)
  )

if (anyNA(analysis_data$treatment) || !all(analysis_data$treatment %in% c(0L, 1L))) {
  stop("Treatment could not be coded as 0/1 from TRTCD.")
}

model_variables <- c("age", "sex", "b_ecogct", "diagtype", "hpv", "dstatus",
                     "prior_hn_treatment")
if (anyNA(analysis_data[, model_variables])) {
  stop("Baseline covariates contain missing values after endpoint construction.")
}

if (!all(analysis_data$complete_record_days == analysis_data$cumulative_severe_ae_duration)) {
  stop("Complete date intervals do not reproduce the prior record-duration endpoint.")
}

arm_quality_summary <- analysis_data %>%
  group_by(treatment_label) %>%
  summarise(
    Patients = n(),
    `Mean non-overlapping severe-AE days` = mean(nonoverlap_severe_ae_days),
    `Median non-overlapping severe-AE days` = median(nonoverlap_severe_ae_days),
    `Zero non-overlapping severe-AE days` = sum(nonoverlap_severe_ae_days == 0),
    `Mean record-days` = mean(complete_record_days),
    `Mean overlap days removed` = mean(overlap_days_removed),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = -treatment_label, names_to = "Metric", values_to = "Value"
  ) %>%
  mutate(Metric = paste(treatment_label, Metric, sep = ": ")) %>%
  select(Metric, Value)

quality_summary <- bind_rows(
  tibble(
    Metric = c(
      "Randomized patients", "Patients with at least one AE record",
      "On-study grade 3+ AE records", "Patients with on-study grade 3+ AE",
      "Grade 3+ records with complete raw start/end days",
      "Grade 3+ records with complete imputed start/end days",
      "Grade 3+ records with unresolved imputed start or end day",
      "Unresolved records missing the end day only",
      "Patients with at least one unresolved record",
      "Patients with no complete qualifying interval, assigned zero",
      "Patients with zero non-overlapping severe-AE days",
      "Patients with at least one day removed through interval overlap",
      "Total record-days before interval union",
      "Total distinct calendar days after interval union",
      "Total overlap days removed"
    ),
    Value = c(
      nrow(analysis_data), n_distinct(ae$SUBJID), nrow(qualifying_events),
      n_distinct(qualifying_events$SUBJID),
      sum(is.finite(qualifying_events$AESTDY) & is.finite(qualifying_events$AEENDY)),
      nrow(complete_events),
      sum(!qualifying_events$complete_imputed_interval),
      sum(with(qualifying_events[!qualifying_events$complete_imputed_interval, ],
               is.finite(AESTDYI) & !is.finite(AEENDYI))),
      n_distinct(qualifying_events$SUBJID[!qualifying_events$complete_imputed_interval]),
      sum(patient_completeness$complete_records == 0L),
      sum(analysis_data$nonoverlap_severe_ae_days == 0),
      sum(analysis_data$overlap_days_removed > 0),
      sum(analysis_data$complete_record_days),
      sum(analysis_data$nonoverlap_severe_ae_days),
      sum(analysis_data$overlap_days_removed)
    )
  ),
  arm_quality_summary
)

write_csv(analysis_data, file.path(out_dir, "nonoverlap_calendar_day_analysis_data.csv"))
write_csv(quality_summary, file.path(out_dir, "nonoverlap_calendar_day_data_quality.csv"))
write_csv(
  analysis_data %>%
    select(SUBJID, treatment_label, cumulative_severe_ae_duration,
           complete_record_days, nonoverlap_severe_ae_days,
           overlap_days_removed, complete_severe_ae_records),
  file.path(out_dir, "nonoverlap_calendar_day_patient_overlap.csv")
)

plot_data <- analysis_data %>%
  select(treatment_label, cumulative_severe_ae_duration, nonoverlap_severe_ae_days) %>%
  tidyr::pivot_longer(
    cols = c(cumulative_severe_ae_duration, nonoverlap_severe_ae_days),
    names_to = "Endpoint", values_to = "Days"
  ) %>%
  mutate(
    Endpoint = recode(
      Endpoint,
      cumulative_severe_ae_duration = "Summed AE-record days",
      nonoverlap_severe_ae_days = "Distinct calendar days"
    )
  )

png(file.path(out_dir, "endpoint_distribution_by_treatment.png"), width = 1000, height = 650)
print(
  ggplot(plot_data, aes(x = treatment_label, y = Days, fill = treatment_label)) +
    geom_boxplot(outlier.alpha = 0.35) +
    facet_wrap(~ Endpoint, scales = "free_y") +
    scale_y_continuous(trans = "sqrt") +
    labs(
      x = NULL, y = "Days (square-root scale)", fill = NULL,
      title = "Severe adverse-event burden by endpoint and treatment arm"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
)
dev.off()

cat("Endpoint construction completed.\n")
cat("Patients:", nrow(analysis_data), "\n")
cat("Complete severe-AE intervals:", nrow(complete_events), "of", nrow(qualifying_events), "\n")
cat("Distinct calendar-day total:", sum(analysis_data$nonoverlap_severe_ae_days), "\n")
