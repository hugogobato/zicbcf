# Shared definitions for the MEPS 2023 dental-expenditure applied study.
#
# Every fit in this directory draws its cleaning rules, design matrix, subgroup
# definitions and survey-design handling from this file. Three defects in the
# first version of the analysis are corrected here.
#
# 1. CATEGORICAL ENCODING. The earlier scripts built the design matrix with
#    as.matrix(df[, covars]), a direct numeric coercion. Race, poverty category,
#    census region and marital status therefore entered the forests as single
#    continuous columns whose values are arbitrary category codes, so a split at
#    REGION23 <= 2 separated Northeast and Midwest from South and West purely
#    because of the order in which the codes were assigned, and a split on race
#    ordered the categories along the code axis. All nominal variables are now
#    expanded with model.matrix, matching the convention already used by the
#    oncology scripts.
#
# 2. HISPANIC ETHNICITY. The earlier scripts derived race from RACEV2X, which
#    carries no ethnicity dimension. In this sample that placed 3,758 Hispanic
#    respondents inside the "White" category and 144 inside "Black". Race and
#    ethnicity are now taken from RACETHX, the standard MEPS analytic variable,
#    whose first level is Hispanic.
#
# 3. SURVEY DESIGN. MEPS is a stratified, clustered, unequal-probability
#    sample. The earlier scripts used no weights, no strata and no primary
#    sampling units, and cut income tertiles at unweighted quantiles, so the
#    reported effect was a sample average over an oversampled design rather than
#    a United States population quantity, and its interval ignored both
#    stratification and clustering. PERWT23F, VARSTR and VARPSU are used here.

MEPS_COVARIATES <- c("AGE23X", "SEX", "RACETHX", "FAMINC23", "POVCAT23",
                     "REGION23", "MARRY23X")
MEPS_DESIGN_VARIABLES <- c("PERWT23F", "VARSTR", "VARPSU")
MEPS_MIN_SUBGROUP_N <- 30L

meps_load <- function(path = "applied_study/h251.csv") {
  if (!file.exists(path)) {
    stop("Missing MEPS file: ", path, ". Run applied_study/download_meps_h251.R first.")
  }
  keep <- c(MEPS_COVARIATES, MEPS_DESIGN_VARIABLES, "DNTINS23_M23", "DVTEXP23")
  df <- read.csv(path)
  missing_columns <- setdiff(keep, names(df))
  if (length(missing_columns) > 0L) {
    stop("MEPS file is missing: ", paste(missing_columns, collapse = ", "))
  }
  df <- df[, keep]

  # Negative codes are MEPS out-of-scope, refused, don't-know and not-ascertained
  # values, not measurements.
  df <- df[df$DNTINS23_M23 > 0, ]
  for (column in MEPS_COVARIATES) df <- df[df[[column]] >= 0, ]
  df <- df[df$DVTEXP23 >= 0, ]

  df$z <- as.integer(df$DNTINS23_M23 == 1)
  df$y <- df$DVTEXP23
  rownames(df) <- NULL
  df
}

meps_design_matrix <- function(df) {
  X <- model.matrix(
    ~ AGE23X + FAMINC23 + factor(SEX) + factor(RACETHX) + factor(POVCAT23) +
      factor(REGION23) + factor(MARRY23X),
    data = df
  )[, -1, drop = FALSE]
  storage.mode(X) <- "double"
  colnames(X) <- gsub("factor\\(|\\)", "", colnames(X))
  if (anyNA(X)) stop("Design matrix contains missing values.")
  if (qr(X)$rank < ncol(X)) stop("Design matrix is rank deficient.")
  X
}

# Weighted quantiles, used so that the income tertile cutpoints describe the
# population rather than the sample.
meps_weighted_quantile <- function(x, weights, probs) {
  keep <- weights > 0
  x <- x[keep]
  weights <- weights[keep]
  order_index <- order(x)
  x <- x[order_index]
  weights <- weights[order_index]
  cumulative <- cumsum(weights) / sum(weights)
  vapply(probs, function(p) {
    if (p <= 0) return(x[1L])
    if (p >= 1) return(x[length(x)])
    x[which(cumulative >= p)[1L]]
  }, numeric(1))
}

meps_subgroup_labels <- function(df) {
  income_breaks <- unique(c(
    min(df$FAMINC23),
    meps_weighted_quantile(df$FAMINC23, df$PERWT23F, c(1 / 3, 2 / 3)),
    max(df$FAMINC23)
  ))
  list(
    Sex = factor(df$SEX, levels = c(1, 2), labels = c("Male", "Female")),
    `Race and ethnicity` = factor(
      df$RACETHX, levels = 1:5,
      labels = c("Hispanic", "White, non-Hispanic", "Black, non-Hispanic",
                 "Asian, non-Hispanic", "Other or multiple, non-Hispanic")),
    `Age group` = cut(df$AGE23X, breaks = c(-Inf, 17, 34, 49, 64, Inf),
                      labels = c("0-17", "18-34", "35-49", "50-64", "65+")),
    `Poverty cat.` = factor(df$POVCAT23, levels = 1:5,
                            labels = c("Poor", "Near-poor", "Low-income",
                                       "Middle-income", "High-income")),
    `Family income (weighted tertile)` = cut(
      df$FAMINC23, breaks = income_breaks, include.lowest = TRUE,
      labels = c("Low-income tertile", "Middle-income tertile",
                 "High-income tertile")[seq_len(length(income_breaks) - 1L)]),
    Region = factor(df$REGION23, levels = 1:4,
                    labels = c("Northeast", "Midwest", "South", "West")),
    `Marital status` = factor(
      ifelse(df$MARRY23X == 1, "Married", "Not married"),
      levels = c("Married", "Not married")),
    `Marital status (detailed)` = factor(
      df$MARRY23X, levels = 1:6,
      labels = c("Married", "Widowed", "Divorced", "Separated",
                 "Never married", "Under 16"))
  )
}

# ---------------------------------------------------------------------------
# Design-based variance of a weighted (domain) mean
#
# Taylor linearization for a stratified multistage design. The linearized value
# for unit i is u_i = w_i d_i (tau_i - theta_d) / sum_j w_j d_j, where d_i is the
# domain indicator; units outside the domain contribute zero but still belong to
# their primary sampling unit. Variance is then the usual between-PSU
# within-stratum sum of squares. Every stratum in the 2023 file carries at least
# two primary sampling units, so no stratum needs collapsing.
# ---------------------------------------------------------------------------
meps_design_variance <- function(values, weights, strata, psu, domain = NULL) {
  if (is.null(domain)) domain <- rep(TRUE, length(values))
  effective <- weights * as.numeric(domain)
  total <- sum(effective)
  if (total <= 0) return(NA_real_)
  theta <- sum(effective * values) / total
  linearized <- effective * (values - theta) / total

  key <- paste(strata, psu, sep = "|")
  psu_total <- tapply(linearized, key, sum)
  psu_stratum <- vapply(strsplit(names(psu_total), "|", fixed = TRUE),
                        function(parts) parts[1L], character(1))
  variance <- 0
  for (stratum in unique(psu_stratum)) {
    totals <- psu_total[psu_stratum == stratum]
    n_psu <- length(totals)
    if (n_psu < 2L) next
    variance <- variance + n_psu / (n_psu - 1L) * sum((totals - mean(totals))^2)
  }
  variance
}

# Design-based variance of the difference between two domain means. The
# linearized variable is the difference of the two domains' linearized values,
# so the between-domain covariance induced by sharing primary sampling units is
# retained rather than assumed away.
meps_design_variance_contrast <- function(values, weights, strata, psu,
                                          domain_a, domain_b) {
  linearize <- function(domain) {
    effective <- weights * as.numeric(domain)
    total <- sum(effective)
    if (total <= 0) return(rep(NA_real_, length(values)))
    theta <- sum(effective * values) / total
    effective * (values - theta) / total
  }
  linearized <- linearize(domain_a) - linearize(domain_b)
  if (anyNA(linearized)) return(NA_real_)

  key <- paste(strata, psu, sep = "|")
  psu_total <- tapply(linearized, key, sum)
  psu_stratum <- vapply(strsplit(names(psu_total), "|", fixed = TRUE),
                        function(parts) parts[1L], character(1))
  variance <- 0
  for (stratum in unique(psu_stratum)) {
    totals <- psu_total[psu_stratum == stratum]
    n_psu <- length(totals)
    if (n_psu < 2L) next
    variance <- variance + n_psu / (n_psu - 1L) * sum((totals - mean(totals))^2)
  }
  variance
}

# Posterior summary of a draw vector, with an optional design-based variance
# component added to the posterior variance.
#
# The reported interval is the posterior mean plus or minus 1.96 times the
# square root of the sum of the two variances. This treats sampling variability
# and posterior uncertainty as approximately independent, which is the standard
# decomposition when a model estimate is carried onto a complex-survey design.
# Clustering usually widens the interval and stratification usually narrows it,
# so the direction of the net change is not knowable in advance.
meps_posterior_summary <- function(draws, design_variance = NA_real_) {
  design_variance <- unname(design_variance)
  quantiles <- quantile(draws, c(0.025, 0.975))
  posterior_variance <- var(draws)
  out <- c(
    Estimate = mean(draws),
    CI_low = unname(quantiles[1L]),
    CI_high = unname(quantiles[2L]),
    P_gt_0 = mean(draws > 0),
    Posterior_SD = sqrt(posterior_variance),
    Design_SE = sqrt(design_variance),
    Total_SE = sqrt(posterior_variance + design_variance)
  )
  out["CI_low_design"] <- out["Estimate"] - 1.96 * out["Total_SE"]
  out["CI_high_design"] <- out["Estimate"] + 1.96 * out["Total_SE"]
  out
}

# Weighted average of unit-level effects within each posterior draw.
meps_weighted_draw_means <- function(cate_draws, weights, index = NULL) {
  if (is.null(index)) index <- seq_len(ncol(cate_draws))
  w <- weights[index]
  if (sum(w) <= 0) return(rep(NA_real_, nrow(cate_draws)))
  as.numeric(cate_draws[, index, drop = FALSE] %*% w) / sum(w)
}
