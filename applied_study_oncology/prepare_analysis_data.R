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
#                       It enters the adjustment set (see oncology_common.R).
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
  library(readr)
})

outdir <- "applied_study_oncology"
patient_file <- file.path(outdir, "headneck_265_severe_ae_patient_level.csv")
out_file <- file.path(outdir, "zic_bcf_headneck_analysis_data.csv")
if (!file.exists(patient_file)) stop("Missing patient-level file: ", patient_file)

patients <- read_csv(patient_file, show_col_types = FALSE)

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
    treatment_duration = TRTDUR
  )

if (anyNA(analysis_data$treatment)) stop("Treatment could not be coded from TRTCD.")
model_variables <- c("age", "b_ecogct", "sex", "diagtype", "hpv", "dstatus",
                     "prior_hn_treatment")
if (anyNA(analysis_data[, model_variables])) {
  stop("Adjustment-set covariates contain missing values.")
}
if (any(analysis_data$treatment_duration <= 0)) {
  stop("Treatment duration must be strictly positive.")
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
cat("\nTreatment duration (TRTDUR) by arm, days:\n")
print(analysis_data %>% group_by(treatment_label) %>%
        summarise(n = n(), mean = mean(treatment_duration),
                  median = median(treatment_duration),
                  max = max(treatment_duration), .groups = "drop"))
