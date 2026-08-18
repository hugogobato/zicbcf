## ---------------------------------------------------------------------------
## Provenance audit of the SPECTRUM analytic cohort: what the supplied extract
## can and cannot establish about the step from the parent trial's 657
## randomized patients to the 520 subjects distributed by Project Data Sphere.
##
## This script reads the raw core-variable file as delivered
## (corevar.sas7bdat) rather than any derived analysis file, because the whole
## point is to characterise the release itself. It answers three questions.
##
##   1. Where does the reduction happen? Not in any code in this repository:
##      the delivered file already contains exactly 520 subject records. The
##      release therefore applied a rule we cannot see, and the extract carries
##      no disposition or exclusion-reason field with which to recover it.
##
##   2. What did the release condition on? The Y/N flag EVALSAFE ("Included in
##      Safety Set?", per the data definition table) and the actual-treatment
##      field ATRT are checked for constancy. If every delivered subject is
##      EVALSAFE = Y with a non-missing ATRT, the release is confined to dosed
##      patients, which accounts for part of the reduction but not for all of
##      it.
##
##   3. Is what remains arm-differential? This is the question that matters for
##      the causal reading, because randomization attaches to the 657 and not
##      automatically to a subset of them. Three tests are run: per-covariate
##      arm balance, an omnibus likelihood-ratio test of assignment on the full
##      baseline covariate set, and the same test restricted to the three IVRS
##      randomization stratification variables. Balance is necessary and not
##      sufficient: a rule operating on a post-randomization variable would not
##      show up in any of these, and the write-up says so.
##
## A fourth, weaker diagnostic reconstructs the subject-numbering scheme.
## SUBJID is "22" concatenated with SITEID and a three-digit within-site
## sequence, so interior gaps in that sequence identify subjects who were
## enrolled at a retained centre and are absent from the release. Subject
## numbers are normally assigned at screening, so the gap count also absorbs
## screen failures who were never randomized; it bounds rather than estimates
## the number of randomized patients removed, and it is reported that way.
##
## Outputs (all withheld from version control, see .gitignore):
##   cohort_provenance_summary.csv
##   cohort_provenance_balance.csv
##   cohort_provenance_site_gaps.csv
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(haven)
})

script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f) == 1L) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())

core_path <- file.path(script_dir, "corevar.sas7bdat")
ae_path <- file.path(script_dir, "New_Data", "ae.sas7bdat")

if (!file.exists(core_path)) {
  stop("corevar.sas7bdat not found. This script requires the raw Project Data ",
       "Sphere release and can only be run by a registered Authorized User.")
}

core <- haven::read_sas(core_path)
core$arm <- ifelse(core$TRT == "Chemotherapy", "Chemotherapy", "Panitumumab")
core$z <- as.integer(core$arm == "Panitumumab")

## --- 1. Where the reduction happens ----------------------------------------

n_rows <- nrow(core)
n_subjects <- length(unique(core$SUBJID))
n_sites <- length(unique(core$SITEID))
arm_counts <- table(core$arm)

## --- 2. What the release conditioned on ------------------------------------

evalsafe_levels <- sort(unique(as.character(core$EVALSAFE)))
n_atrt_missing <- sum(is.na(core$ATRT) | core$ATRT == "")
n_trtdur_missing <- sum(is.na(core$TRTDUR))

## --- 3. Whether the remainder is arm-differential --------------------------

balance_vars <- c(
  AGECAT = "Age category", SEX = "Sex", RACE = "Race",
  B_ECOGCT = "Baseline ECOG category (CRF)",
  BECGICT = "Baseline ECOG category (IVRS stratum)",
  DIAGTYPE = "Primary tumour diagnosis",
  TUMCAT = "Tumour category (CRF)",
  TUMCATI = "Tumour category (IVRS stratum)",
  PRHNTRTC = "Prior SCCHN treatment (CRF composite)",
  PRHNTRTI = "Prior SCCHN treatment (IVRS stratum)",
  DSTATUS = "Disease status", HPV = "HPV status"
)

balance_rows <- lapply(names(balance_vars), function(v) {
  x <- as.character(core[[v]])
  x[is.na(x) | x == ""] <- "(not recorded)"
  tb <- table(x, core$arm)
  p <- tryCatch(suppressWarnings(stats::chisq.test(tb)$p.value),
                error = function(e) NA_real_)
  data.frame(
    variable = v,
    label = unname(balance_vars[v]),
    level = rownames(tb),
    chemotherapy_n = as.integer(tb[, "Chemotherapy"]),
    chemotherapy_pct = round(100 * tb[, "Chemotherapy"] / sum(tb[, "Chemotherapy"]), 1),
    panitumumab_n = as.integer(tb[, "Panitumumab"]),
    panitumumab_pct = round(100 * tb[, "Panitumumab"] / sum(tb[, "Panitumumab"]), 1),
    chisq_p = round(p, 3),
    stringsAsFactors = FALSE
  )
})
balance <- do.call(rbind, balance_rows)

age_test <- stats::t.test(AGE ~ arm, data = core)
balance <- rbind(balance, data.frame(
  variable = "AGE", label = "Age at screening (years, mean)", level = "(continuous)",
  chemotherapy_n = sum(core$arm == "Chemotherapy"),
  chemotherapy_pct = round(mean(core$AGE[core$arm == "Chemotherapy"]), 2),
  panitumumab_n = sum(core$arm == "Panitumumab"),
  panitumumab_pct = round(mean(core$AGE[core$arm == "Panitumumab"]), 2),
  chisq_p = round(age_test$p.value, 3), stringsAsFactors = FALSE
))

min_p <- min(balance$chisq_p, na.rm = TRUE)
## HPV is reported separately because its "not recorded" level reflects whether
## tumour tissue was provided and assayed, which is not a baseline covariate in
## the randomization sense even though the biology it measures is.
min_p_excl_hpv <- min(balance$chisq_p[balance$variable != "HPV"], na.rm = TRUE)

core$RACE2 <- ifelse(core$RACE == "White or Caucasian", "White",
                     ifelse(core$RACE == "Asian", "Asian", "Other or not recorded"))

fit_full <- stats::glm(
  z ~ AGE + factor(SEX) + factor(RACE2) + factor(B_ECOGCT) + factor(BECGICT) +
    factor(DIAGTYPE) + factor(DSTATUS) + factor(PRHNTRTC),
  family = stats::binomial(), data = core
)
fit_null <- stats::glm(z ~ 1, family = stats::binomial(), data = core)
fit_ivrs <- stats::glm(
  z ~ factor(BECGICT) + factor(TUMCATI) + factor(PRHNTRTI),
  family = stats::binomial(), data = core
)

lr_test <- function(fit, null) {
  stat <- as.numeric(2 * (stats::logLik(fit) - stats::logLik(null)))
  df <- fit$rank - null$rank
  c(chisq = stat, df = df, p = stats::pchisq(stat, df, lower.tail = FALSE))
}
lr_full <- lr_test(fit_full, fit_null)
lr_ivrs <- lr_test(fit_ivrs, fit_null)
ps <- stats::predict(fit_full, type = "response")

## --- 4. Subject-numbering diagnostic ---------------------------------------

site_from_id <- substr(core$SUBJID, 3, 6)
numbering_consistent <- all(site_from_id == core$SITEID)
seq_no <- as.integer(substr(core$SUBJID, 7, 9))

site_gaps <- do.call(rbind, lapply(split(seq_no, core$SITEID), function(s) {
  data.frame(n_released = length(s), max_sequence = max(s),
             interior_gaps = max(s) - length(s), stringsAsFactors = FALSE)
}))
site_gaps <- data.frame(site = rownames(site_gaps), site_gaps,
                        row.names = NULL, stringsAsFactors = FALSE)
site_gaps <- site_gaps[order(-site_gaps$interior_gaps, site_gaps$site), ]

total_gaps <- sum(site_gaps$interior_gaps)
sites_with_gaps <- sum(site_gaps$interior_gaps > 0)

## Country-block prefixes, as a note only: the release carries no country field.
prefixes <- sort(unique(substr(core$SITEID, 1, 2)))

## --- 5. Linkage to the adverse-event file ----------------------------------

n_ae_subjects <- NA_integer_
n_no_ae <- NA_integer_
if (file.exists(ae_path)) {
  ae <- haven::read_sas(ae_path, col_select = c("SUBJID"))
  n_ae_subjects <- length(unique(ae$SUBJID))
  n_no_ae <- length(setdiff(core$SUBJID, unique(ae$SUBJID)))
}

## --- Write-out --------------------------------------------------------------

summary_tbl <- data.frame(
  quantity = c(
    "Randomized in the parent trial (published report)",
    "Subject records in the delivered core-variable file",
    "Distinct subjects in the delivered core-variable file",
    "Chemotherapy arm",
    "Panitumumab plus chemotherapy arm",
    "Distinct investigational sites represented",
    "Distinct levels of EVALSAFE in the release",
    "Subjects with missing actual treatment (ATRT)",
    "Subjects with missing treatment duration (TRTDUR)",
    "Smallest arm-balance p-value across baseline covariates",
    "Smallest arm-balance p-value excluding HPV result availability",
    "Omnibus balance LR chi-square (all baseline covariates)",
    "Omnibus balance LR degrees of freedom",
    "Omnibus balance LR p-value",
    "IVRS-stratum-only balance LR chi-square",
    "IVRS-stratum-only balance LR degrees of freedom",
    "IVRS-stratum-only balance LR p-value",
    "Estimated propensity score, minimum",
    "Estimated propensity score, maximum",
    "Estimated propensity score, standard deviation",
    "SUBJID decomposes as 22 + SITEID + within-site sequence",
    "Sites with at least one interior gap in the sequence",
    "Total interior gaps in the within-site sequence",
    "Distinct site-identifier country blocks represented",
    "Distinct subjects in the adverse-event file",
    "Released subjects with no adverse-event record"
  ),
  value = c(
    657, n_rows, n_subjects,
    as.integer(arm_counts["Chemotherapy"]), as.integer(arm_counts["Panitumumab"]),
    n_sites, length(evalsafe_levels), n_atrt_missing, n_trtdur_missing,
    min_p, min_p_excl_hpv,
    round(lr_full["chisq"], 2), lr_full["df"], round(lr_full["p"], 3),
    round(lr_ivrs["chisq"], 2), lr_ivrs["df"], round(lr_ivrs["p"], 3),
    round(min(ps), 3), round(max(ps), 3), round(stats::sd(ps), 3),
    as.integer(numbering_consistent), sites_with_gaps, total_gaps,
    length(prefixes), n_ae_subjects, n_no_ae
  ),
  note = c(
    "Vermorken et al. (2013); not reproducible from the supplied files",
    "corevar.sas7bdat exactly as delivered; the reduction precedes any code here",
    "One row per subject",
    "", "",
    "SITEID values are de-identified",
    paste0("Observed level(s): ", paste(evalsafe_levels, collapse = ", ")),
    "Zero implies every released subject was dosed",
    "Zero implies every released subject was dosed",
    "Chi-square, or t-test for age; larger is more consistent with balance",
    "HPV availability depends on tissue provision, which is not a baseline covariate",
    "Logistic regression of arm on the full baseline covariate set", "", "",
    "Restricted to the three IVRS randomization stratification variables", "", "",
    "Centred near 0.5 under intact randomization", "", "",
    "1 = yes; enables the gap diagnostic below",
    paste0("Of ", nrow(site_gaps), " sites"),
    "Upper bound: subject numbers are assigned at screening, so screen failures also consume numbers",
    paste0("Blocks present: ", paste(prefixes, collapse = ", ")),
    "New_Data/ae.sas7bdat",
    "Each contributes an outcome of zero at every horizon"
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(summary_tbl, file.path(script_dir, "cohort_provenance_summary.csv"), row.names = FALSE)
utils::write.csv(balance, file.path(script_dir, "cohort_provenance_balance.csv"), row.names = FALSE)
utils::write.csv(site_gaps, file.path(script_dir, "cohort_provenance_site_gaps.csv"), row.names = FALSE)

cat("\n=== Cohort provenance audit ===\n")
print(summary_tbl[, c("quantity", "value")], right = FALSE)
cat("\nBalance table written to cohort_provenance_balance.csv (",
    nrow(balance), " rows).\n", sep = "")
cat("Site-level sequence gaps written to cohort_provenance_site_gaps.csv.\n")
