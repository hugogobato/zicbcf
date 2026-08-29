# Reconcile the second Project Data Sphere download with the AE-derived
# oncology analysis and quantify what the new domains do and do not identify.
#
# The script is descriptive. It does not replace the randomized estimand with
# a rate, condition on treatment duration, or claim that last contact is a
# complete safety-assessment window. It uses LASTCTDY as an observed upper
# bound for a censoring diagnostic and DSDY as the recorded safety-follow-up
# day, retaining their distinct meanings.

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(readr)
  library(tidyr)
})

script_file <- sub("^--file=", "",
                   commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))])
if (length(script_file) != 1L) stop("Run this script with Rscript.")

domain_dir <- dirname(normalizePath(script_file))
study_dir <- dirname(domain_dir)
out_dir <- file.path(domain_dir, "analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_domain <- function(name) {
  read_sas(file.path(domain_dir, name), .name_repair = "minimal") %>%
    mutate(SUBJID = as.character(SUBJID))
}

demo <- read_domain("demo.sas7bdat")
disposit <- read_domain("disposit.sas7bdat")
radiotx <- read_domain("radiotx.sas7bdat")
lesion <- read_domain("lesion.sas7bdat")
ae <- read_sas(file.path(study_dir, "New_Data", "ae.sas7bdat"),
               .name_repair = "minimal") %>%
  mutate(
    SUBJID = as.character(SUBJID),
    qualifying = AESEVCD >= 3 & (is.na(AEPRIOR) | AEPRIOR != "Y"),
    complete_interval = is.finite(AESTDYI) & is.finite(AEENDYI)
  )

patients <- demo %>%
  transmute(
    SUBJID,
    arm = if_else(TRTCD == 2, "Panitumumab + chemotherapy", "Chemotherapy"),
    treatment_code = TRTCD,
    treatment_duration = TRTDUR,
    prior_hn_treatment = PRHNTRTC,
    prior_radiotherapy = PRRADIO
  ) %>%
  left_join(disposit %>% select(SUBJID, DSSTAT, DSDY, LASTOSDY, LASTCTDY,
                                AFUP, PFUP, DSEOS, EOIPDY), by = "SUBJID")

if (nrow(patients) != 520L || n_distinct(patients$SUBJID) != 520L ||
    anyNA(patients$LASTCTDY) || anyNA(patients$DSDY)) {
  stop("The v2 domains do not provide one complete row per analytic subject.")
}

if (!all(table(patients$arm) == 260L)) {
  stop("The v2 demo file is not balanced 260/260 by assigned arm.")
}

write_csv(
  patients %>%
    count(arm, prior_radiotherapy) %>%
    group_by(arm) %>%
    mutate(percent = 100 * n / sum(n)) %>%
    ungroup(),
  file.path(out_dir, "prior_radiotherapy_by_arm.csv")
)

radio_subject <- radiotx %>%
  group_by(SUBJID) %>%
  summarise(
    radiotx_any = any(RAANY == "Y"),
    radiotx_records = n(),
    radiotx_sites = n_distinct(RASITE[RAANY == "Y" & !is.na(RASITE) & RASITE != ""]),
    .groups = "drop"
  )
write_csv(
  patients %>%
    select(SUBJID, arm, prior_radiotherapy) %>%
    left_join(radio_subject, by = "SUBJID") %>%
    mutate(radiotx_any = coalesce(radiotx_any, FALSE)) %>%
    count(prior_radiotherapy, radiotx_any),
  file.path(out_dir, "prior_radiotherapy_concordance.csv")
)

lesion_subject <- lesion %>%
  group_by(SUBJID) %>%
  summarise(
    lesion_prior_rt_yes = any(LSPIRRNY == "Y"),
    lesion_prior_rt_no = any(LSPIRRNY == "N"),
    lesion_records = n(),
    .groups = "drop"
  )
write_csv(
  patients %>%
    select(SUBJID, arm, prior_radiotherapy) %>%
    left_join(lesion_subject, by = "SUBJID") %>%
    mutate(across(c(lesion_prior_rt_yes, lesion_prior_rt_no), ~ coalesce(.x, FALSE))) %>%
    count(prior_radiotherapy, lesion_prior_rt_yes, lesion_prior_rt_no),
  file.path(out_dir, "lesion_prior_radiotherapy_summary.csv")
)

write_csv(
  patients %>%
    group_by(arm) %>%
    summarise(
      n = n(),
      safety_followup_completed = sum(DSSTAT == "Completed safety FUP"),
      safety_followup_completed_percent = 100 * mean(DSSTAT == "Completed safety FUP"),
      mean_safety_followup_day = mean(DSDY),
      median_safety_followup_day = median(DSDY),
      mean_last_contact_day = mean(LASTCTDY),
      median_last_contact_day = median(LASTCTDY),
      last_contact_before_84 = sum(LASTCTDY < 84),
      last_contact_before_120 = sum(LASTCTDY < 120),
      reaches_84_by_last_contact = sum(LASTCTDY >= 84),
      reaches_120_by_last_contact = sum(LASTCTDY >= 120),
      .groups = "drop"
    ),
  file.path(out_dir, "followup_by_arm.csv")
)

qualifying <- ae %>% filter(qualifying)
complete <- qualifying %>% filter(complete_interval)
complete_bounds <- complete %>%
  left_join(patients %>% select(SUBJID, DSDY, LASTCTDY), by = "SUBJID")
incomplete_bounds <- qualifying %>%
  filter(!complete_interval) %>%
  left_join(patients %>% select(SUBJID, DSDY, LASTCTDY), by = "SUBJID")
write_csv(
  tibble(
    quantity = c(
      "Qualifying on-study grade 3+ records",
      "Complete imputed intervals",
      "Unresolved intervals",
      "Complete intervals ending after DSDY",
      "Complete intervals ending after LASTCTDY",
      "Complete intervals starting after LASTCTDY",
      "Unresolved intervals with a start after DSDY",
      "Unresolved intervals with a start after LASTCTDY"
    ),
    value = c(
      nrow(qualifying), nrow(complete), sum(!qualifying$complete_interval),
      sum(complete_bounds$AEENDYI > complete_bounds$DSDY),
      sum(complete_bounds$AEENDYI > complete_bounds$LASTCTDY),
      sum(complete_bounds$AESTDYI > complete_bounds$LASTCTDY),
      sum(incomplete_bounds$AESTDYI > incomplete_bounds$DSDY, na.rm = TRUE),
      sum(incomplete_bounds$AESTDYI > incomplete_bounds$LASTCTDY, na.rm = TRUE)
    )
  ),
  file.path(out_dir, "followup_event_boundary_diagnostics.csv")
)

union_days <- function(start_day, end_day, horizon = Inf) {
  start_day <- pmax(start_day, 1)
  end_day <- pmin(end_day, horizon)
  keep <- start_day <= end_day
  if (!any(keep)) return(0)
  start_day <- start_day[keep]
  end_day <- end_day[keep]
  order_index <- order(start_day, end_day)
  start_day <- start_day[order_index]
  end_day <- end_day[order_index]
  current_start <- start_day[1L]
  current_end <- end_day[1L]
  total <- 0
  if (length(start_day) > 1L) {
    for (i in seq.int(2L, length(start_day))) {
      if (start_day[i] <= current_end + 1L) {
        current_end <- max(current_end, end_day[i])
      } else {
        total <- total + current_end - current_start + 1L
        current_start <- start_day[i]
        current_end <- end_day[i]
      }
    }
  }
  total + current_end - current_start + 1L
}

complete_with_followup <- complete %>%
  left_join(patients %>% select(SUBJID, arm, LASTCTDY), by = "SUBJID")
horizons <- c(84, 90, 120)
# Last-contact capping applies only to intervals with a complete boundary, so
# the horizon diagnostics below inherit, in bounded form, the forced-zero
# problem created by the unresolved-end-day records documented upstream.
observed_window <- complete_with_followup %>%
  group_by(SUBJID, arm, LASTCTDY) %>%
  summarise(
    last_contact_capped_days = union_days(AESTDYI, AEENDYI, LASTCTDY),
    h84_observed_days = union_days(AESTDYI, AEENDYI, pmin(84, LASTCTDY)),
    h90_observed_days = union_days(AESTDYI, AEENDYI, pmin(90, LASTCTDY)),
    h120_observed_days = union_days(AESTDYI, AEENDYI, pmin(120, LASTCTDY)),
    .groups = "drop"
  )
observed_patient <- patients %>%
  left_join(observed_window, by = c("SUBJID", "arm", "LASTCTDY")) %>%
  mutate(across(c(last_contact_capped_days,
                  all_of(paste0("h", horizons, "_observed_days"))),
                ~ coalesce(.x, 0)))

truncation_rows <- lapply(c(Inf, horizons), function(h) {
  label <- if (is.infinite(h)) "Last-contact-capped" else paste0(h, " days")
  variable <- if (is.infinite(h)) "last_contact_capped_days" else paste0("h", h, "_observed_days")
  observed_patient %>%
    group_by(Arm = arm) %>%
    summarise(
      Horizon = label,
      Mean_observed_days = mean(.data[[variable]]),
      N_reaching_horizon = if (is.infinite(h)) NA_integer_ else sum(LASTCTDY >= h),
      Mean_among_reaching = if (is.infinite(h)) NA_real_ else
        mean(.data[[variable]][LASTCTDY >= h]),
      .groups = "drop"
    )
})
write_csv(bind_rows(truncation_rows), file.path(out_dir, "last_contact_horizon_diagnostic.csv"))

## The per-patient endpoint file is a derivative with one row per subject and is
## therefore never written inside the versioned repository: it goes to a secure
## directory outside it (override with the ZICBCF_SECURE_DIR environment
## variable). Only aggregate outputs are written to out_dir above.
secure_dir <- Sys.getenv("ZICBCF_SECURE_DIR",
                         unset = file.path(path.expand("~"), "secure_zicbcf_outputs",
                                           "New_Data_v2"))
dir.create(secure_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(
  observed_patient %>% select(SUBJID, arm, LASTCTDY,
                              last_contact_capped_days,
                              all_of(paste0("h", horizons, "_observed_days"))),
  file.path(secure_dir, "last_contact_capped_patient_endpoint.csv")
)
cat("Wrote v2 domain diagnostics to", out_dir, "\n")
cat("Wrote per-patient endpoint file outside the repository to", secure_dir, "\n")
