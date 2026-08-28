# Shared definitions for the SPECTRUM oncology applied study.
#
# Every fit in this directory (ZIC-BCF-Smear and the Gamma-hurdle benchmark, on
# both the summed-record-duration and the distinct-calendar-day endpoint) draws
# its adjustment set and its subgroup definitions from this file, so that the
# specifications cannot silently drift apart.
#
# Two changes relative to the first version of the analysis are encoded here.
#
# 1. Tumor category is dropped from the adjustment set. TUMCAT is a
#    deterministic recode of DIAGTYPE in this extract (larynx and oropharynx map
#    to Y; hypopharynx and oral cavity map to N), so including both puts an
#    exact linear combination of existing columns into the design matrix.
#
# 2. The now-available patient-level prior radiotherapy indicator (PRRADIO) is
#    used in the adjustment set. Prior therapy to the head and neck is a
#    principal determinant of mucosal and cutaneous toxicity, which is the
#    outcome being modeled. PRHNTRTC is retained in the analysis data as the
#    broader composite but is not entered alongside PRRADIO, in the adjustment
#    set or the subgroup lists: the two indicators differ for only four
#    patients, so carrying both adds a near-duplicate contrast to every
#    subgroup grid in an analysis that is already unadjusted for multiplicity,
#    and their near-redundancy destabilized the Gamma benchmark.

ONC_COVARIATE_FORMULA <- ~ age + b_ecogct + sex + diagtype + hpv + dstatus +
  prior_radiotherapy

# Minimum subgroup size for a reported subgroup effect or pairwise contrast.
ONC_MIN_SUBGROUP_N <- 30L

onc_design_matrix <- function(df) {
  required <- c("age", "b_ecogct", "sex", "diagtype", "hpv", "dstatus",
                "prior_radiotherapy")
  missing_columns <- setdiff(required, names(df))
  if (length(missing_columns) > 0L) {
    stop("Analysis data is missing: ", paste(missing_columns, collapse = ", "),
         ". Run prepare_analysis_data.R first.")
  }
  X <- model.matrix(ONC_COVARIATE_FORMULA, data = df)[, -1, drop = FALSE]
  storage.mode(X) <- "double"
  if (anyNA(X)) stop("Design matrix contains missing values.")
  if (qr(X)$rank < ncol(X)) {
    stop("Design matrix is rank deficient; a covariate is an exact linear ",
         "combination of the others.")
  }
  X
}

# HPV status is reported with the absent category named "Not recorded" rather
# than "Unknown". The category is created by coalescing missing baseline HPV
# results, so any contrast against it is a missingness contrast and not a
# comparison of two measured biological states.
onc_subgroup_labels <- function(df) {
  list(
    Sex = factor(df$sex, levels = c("Male", "Female")),
    `ECOG status` = factor(df$b_ecogct, levels = c(0, 1),
                           labels = c("ECOG 0", "ECOG 1")),
    `Diagnosis site` = factor(df$diagtype),
    `HPV status` = factor(df$hpv, levels = c("Negative", "Positive", "Unknown"),
                          labels = c("Negative", "Positive", "Not recorded")),
    `Disease status` = factor(df$dstatus,
                              levels = c("newly diagnosed", "recurrent")),
    # PRHNTRTC ("Prior SCCHN treatment") is deliberately absent here: it
    # differs from prior_radiotherapy for four patients, so its contrast grid
    # duplicates the Prior-radiotherapy grid almost exactly.
    `Prior radiotherapy` = factor(df$prior_radiotherapy, levels = c("N", "Y"),
                                  labels = c("No prior radiotherapy",
                                             "Prior radiotherapy")),
    `Age group` = cut(df$age, breaks = c(-Inf, 49, 59, 69, Inf),
                      labels = c("<50", "50-59", "60-69", "70+"))
  )
}

onc_posterior_summary <- function(draws) {
  ci <- quantile(draws, c(0.025, 0.975))
  c(Estimate = mean(draws), CI_low = unname(ci[1L]), CI_high = unname(ci[2L]),
    P_gt_0 = mean(draws > 0))
}

onc_subgroup_table <- function(cate_draws, labels, model, estimand) {
  rows <- list(data.frame(
    Model = model, Estimand = estimand, Covariate = "Overall", Level = "All",
    N = ncol(cate_draws), t(onc_posterior_summary(rowMeans(cate_draws))),
    check.names = FALSE
  ))
  for (group_name in names(labels)) {
    label <- labels[[group_name]]
    for (level in levels(label)) {
      index <- which(label == level)
      if (length(index) < ONC_MIN_SUBGROUP_N) next
      rows[[length(rows) + 1L]] <- data.frame(
        Model = model, Estimand = estimand, Covariate = group_name,
        Level = level, N = length(index),
        t(onc_posterior_summary(rowMeans(cate_draws[, index, drop = FALSE]))),
        check.names = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

onc_contrast_table <- function(cate_draws, labels, model, estimand) {
  rows <- list()
  for (group_name in names(labels)) {
    label <- labels[[group_name]]
    eligible <- levels(label)[table(label) >= ONC_MIN_SUBGROUP_N]
    if (length(eligible) < 2L) next
    comparisons <- combn(eligible, 2L)
    for (column in seq_len(ncol(comparisons))) {
      first <- comparisons[1L, column]
      second <- comparisons[2L, column]
      contrast <- rowMeans(cate_draws[, which(label == first), drop = FALSE]) -
        rowMeans(cate_draws[, which(label == second), drop = FALSE])
      rows[[length(rows) + 1L]] <- data.frame(
        Model = model, Estimand = estimand, Covariate = group_name,
        Contrast = paste(first, "-", second),
        t(onc_posterior_summary(contrast)), check.names = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
