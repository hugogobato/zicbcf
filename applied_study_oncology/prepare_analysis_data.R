# Rebuild the de-identified analytic file used by every model fit in this
# directory from the patient-level screening output.
#
# The only substantive change relative to the previous version of
# zic_bcf_headneck_analysis_data.csv is that two measured patient-level fields
# are now carried through instead of being discarded at this step:
#
#   prior_hn_treatment  PRHNTRTC, prior treatment for squamous-cell carcinoma
#                       of the head and neck. Per the data definition table this
#                       is set to Y if the subject had prior radiotherapy as per
#                       RAHX, or other prior treatment as per TXHX or OTHEPROC.
#                       It is retained for descriptive and subgroup analyses.
#
#   prior_radiotherapy  PRRADIO, the patient-level prior-radiotherapy indicator
#                       supplied in the v2 demo domain. It is the clinically
#                       specific baseline adjustment variable.
#
#   treatment_duration  TRTDUR, defined in the data definition table as
#                       LDOSDT - FDOSDT + 1. This is time on treatment. It is a
#                       post-randomization consequence of assignment and is
#                       therefore NEVER placed in any adjustment set. It is
#                       carried only so that the descriptive exposure analyses
#                       can be run from the analytic file.
#
# The script also writes a cohort-derivation table, because the step from the
# parent trial's randomized population to this analytic extract is otherwise
# invisible to a reader.

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(readr)
})

outdir <- "applied_study_oncology"
patient_file <- file.path(outdir, "headneck_265_severe_ae_patient_level.csv")
out_file <- file.path(outdir, "zic_bcf_headneck_analysis_data.csv")
if (!file.exists(patient_file)) stop("Missing patient-level file: ", patient_file)

patients <- read_csv(patient_file, show_col_types = FALSE) %>%
  mutate(SUBJID = as.character(SUBJID))

# The v2 Project Data Sphere download supplies the patient-level radiotherapy
# indicator and disposition domains that were absent from the first extract.
# Keep the original AE-derived patient file as the outcome source, but verify
# that the new domains link one-to-one to the same 520 subjects.
demo_file <- file.path(outdir, "New_Data_v2", "demo.sas7bdat")
disposit_file <- file.path(outdir, "New_Data_v2", "disposit.sas7bdat")
if (!file.exists(demo_file) || !file.exists(disposit_file)) {
  stop("Missing New_Data_v2/demo.sas7bdat or New_Data_v2/disposit.sas7bdat.")
}
demo_v2 <- read_sas(demo_file, .name_repair = "minimal") %>%
  mutate(SUBJID = as.character(SUBJID)) %>%
  select(SUBJID, PRRADIO)
disposit_v2 <- read_sas(disposit_file, .name_repair = "minimal") %>%
  mutate(SUBJID = as.character(SUBJID)) %>%
  select(SUBJID, DSSTAT, DSDY, LASTOSDY, LASTCTDY, AFUP, PFUP)

if (nrow(demo_v2) != dplyr::n_distinct(demo_v2$SUBJID) ||
    nrow(disposit_v2) != dplyr::n_distinct(disposit_v2$SUBJID) ||
    !setequal(demo_v2$SUBJID, patients$SUBJID) ||
    !setequal(disposit_v2$SUBJID, patients$SUBJID)) {
  stop("New_Data_v2 domains do not link one-to-one to the analytic cohort.")
}

patients <- patients %>%
  left_join(demo_v2, by = "SUBJID") %>%
  left_join(disposit_v2, by = "SUBJID")

required <- c("SUBJID", "AGE", "SEX", "TRTCD", "B_ECOGCT", "DIAGTYPE", "TUMCAT",
              "HPV", "DSTATUS", "PRHNTRTC", "PRHNTRTI", "TRTDUR",
              "severe_ae_total_duration")
missing_columns <- setdiff(required, names(patients))
if (length(missing_columns) > 0L) {
  stop("Patient-level file is missing: ", paste(missing_columns, collapse = ", "))
}
if (nrow(patients) != dplyr::n_distinct(patients$SUBJID)) {
  stop("The patient-level file must contain exactly one row per SUBJID.")
}

analysis_data <- patients %>%
  arrange(SUBJID) %>%
  transmute(
    generated_patient_id = sprintf("PDS_%04d", row_number()),
    treatment = case_when(TRTCD == 2 ~ 1L, TRTCD == 50 ~ 0L, TRUE ~ NA_integer_),
    treatment_label = if_else(TRTCD == 2, "Panitumumab + chemotherapy",
                              "Chemotherapy"),
    cumulative_severe_ae_duration = severe_ae_total_duration,
    any_positive_severe_ae_duration = as.integer(severe_ae_total_duration > 0),
    age = AGE,
    sex = SEX,
    b_ecogct = B_ECOGCT,
    diagtype = DIAGTYPE,
    tumcat = TUMCAT,
    # Absent baseline HPV results are coalesced into a single category. Any
    # contrast against it is a missingness contrast, not a biological one.
    hpv = coalesce(HPV, "Unknown"),
    dstatus = DSTATUS,
    prior_hn_treatment = PRHNTRTC,
    prior_hn_treatment_ivrs = PRHNTRTI,
    prior_radiotherapy = PRRADIO,
    treatment_duration = TRTDUR,
    safety_followup_status = DSSTAT,
    safety_followup_day = DSDY,
    last_on_study_day = LASTOSDY,
    last_contact_day = LASTCTDY,
    actual_followup_weeks = AFUP,
    potential_followup_weeks = PFUP
  )

if (anyNA(analysis_data$treatment)) stop("Treatment could not be coded from TRTCD.")
model_variables <- c("age", "b_ecogct", "sex", "diagtype", "hpv", "dstatus",
                     "prior_hn_treatment", "prior_radiotherapy")
if (anyNA(analysis_data[, model_variables])) {
  stop("Adjustment-set covariates contain missing values.")
}
if (any(analysis_data$treatment_duration <= 0)) {
  stop("Treatment duration must be strictly positive.")
}
if (anyNA(analysis_data$last_contact_day) || any(analysis_data$last_contact_day < 1)) {
  stop("Last-contact study day must be observed and positive.")
}

# Tumor category is retained in the file for descriptive use only. Verify that
# it is the deterministic recode of diagnosis site that the supplement claims,
# so that the justification for excluding it from the design matrix is checked
# rather than asserted.
recode_check <- analysis_data %>%
  count(diagtype, tumcat) %>%
  group_by(diagtype) %>%
  summarise(levels_per_site = n(), .groups = "drop")
if (any(recode_check$levels_per_site != 1L)) {
  stop("Tumor category is not a deterministic function of diagnosis site; ",
       "revisit the exclusion in oncology_common.R.")
}

write_csv(analysis_data, out_file)

cohort_derivation <- tibble(
  Step = c(
    "Randomized in the parent SPECTRUM trial (published report)",
    "Patients in the supplied Project Data Sphere extract, per randomized arm",
    "Patients in the analytic cohort",
    "Analytic cohort: chemotherapy arm",
    "Analytic cohort: panitumumab plus chemotherapy arm",
    "Analytic cohort patients with at least one adverse-event record",
    "Analytic cohort patients with no adverse-event record"
  ),
  Value = c(
    657,
    NA_real_,
    nrow(analysis_data),
    sum(analysis_data$treatment == 0L),
    sum(analysis_data$treatment == 1L),
    NA_real_,
    NA_real_
  ),
  Note = c(
    "Vermorken et al. (2013); not reproducible from the supplied files",
    "Determined by the contents of corevar.sas7bdat in the supplied extract",
    "One row per subject in the supplied core-variable file with a codeable arm",
    "", "",
    "Filled in by New_Data/run_exposure_and_missingness_diagnostics.R",
    "Filled in by New_Data/run_exposure_and_missingness_diagnostics.R"
  )
)
write_csv(cohort_derivation, file.path(outdir, "cohort_derivation.csv"))

cat("Wrote", out_file, "with", nrow(analysis_data), "patients.\n")
cat("Prior SCCHN treatment (PRHNTRTC) by arm:\n")
print(with(analysis_data, table(treatment_label, prior_hn_treatment)))
cat("Prior radiotherapy (PRRADIO) by arm:\n")
print(with(analysis_data, table(treatment_label, prior_radiotherapy)))
cat("\nTreatment duration (TRTDUR) by arm, days:\n")
print(analysis_data %>% group_by(treatment_label) %>%
        summarise(n = n(), mean = mean(treatment_duration),
                  median = median(treatment_duration),
                  max = max(treatment_duration), .groups = "drop"))
