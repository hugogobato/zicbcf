#' Posterior convergence diagnostics for ZIC-BCF-Smear fits
#'
#' Computes standard MCMC convergence diagnostics for a set of monitored scalar
#' functionals of one or more \code{zicbcf_smear()} fits. Three diagnostics are
#' reported: the Gelman--Rubin potential scale reduction factor \eqn{\widehat{R}}
#' (Gelman and Rubin, 1992), the effective sample size, and the
#' Geweke spectral-density test comparing the first and last portions of each
#' chain (Geweke, 1992). All three are obtained from
#' \pkg{coda}.
#'
#' When a single chain is supplied, \eqn{\widehat{R}} is the split-chain variant:
#' each chain is divided into two halves that are then treated as separate
#' sequences, which detects the within-chain drift that a single-chain diagnostic
#' would otherwise miss. When several chains are supplied (recommended: four
#' chains started from different seeds), \eqn{\widehat{R}} is computed across the
#' split halves of all chains.
#'
#' @param fits A single object of class \code{"zicbcf_fit_smear"} or a list of
#'   such objects, one per chain. Chains must have the same number of retained
#'   draws and be fit to the same data.
#' @param n_cate_units Number of unit-level conditional average treatment effects
#'   to monitor in addition to the scalar summaries. Units are selected at evenly
#'   spaced quantiles of the posterior-mean CATE of the first chain, so the
#'   monitored set spans the estimated effect surface rather than an arbitrary
#'   corner of it. Set to zero to monitor only the scalar summaries.
#' @param geweke_frac1,geweke_frac2 Fractions of each chain used for the first
#'   and last windows of the Geweke test.
#'
#' @return A data frame with one row per monitored quantity and columns
#'   \code{quantity}, \code{n_chains}, \code{n_draws}, \code{mean}, \code{rhat},
#'   \code{ess}, \code{ess_per_draw}, \code{mcse}, and \code{geweke_z}
#'   (the maximum absolute Geweke statistic over chains).
#'
#' @references
#' Gelman, A. and Rubin, D. B. (1992). Inference from iterative simulation using
#' multiple sequences. \emph{Statistical Science} 7, 457--472.
#'
#' Geweke, J. (1992). Evaluating the accuracy of sampling-based approaches to
#' the calculation of posterior moments. In \emph{Bayesian Statistics 4}, 169--193.
#'
#' @export
zicbcf_convergence <- function(fits,
                               n_cate_units = 10L,
                               geweke_frac1 = 0.1,
                               geweke_frac2 = 0.5) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for zicbcf_convergence().")
  }
  if (inherits(fits, "zicbcf_fit_smear")) fits <- list(fits)
  if (!is.list(fits) || !length(fits)) {
    stop("'fits' must be a zicbcf_fit_smear object or a non-empty list of them.")
  }
  if (!all(vapply(fits, inherits, logical(1), "zicbcf_fit_smear"))) {
    stop("Every element of 'fits' must be of class 'zicbcf_fit_smear'.")
  }

  n_draws <- vapply(fits, function(f) length(f$ate), integer(1))
  if (length(unique(n_draws)) != 1L) {
    stop("All chains must have the same number of retained draws.")
  }
  n_draws <- n_draws[1L]
  if (n_draws < 8L) stop("Too few retained draws for convergence diagnostics.")

  monitored <- lapply(fits, zicbcf_monitored_scalars)

  # Units are chosen once, from the first chain, so that every chain monitors
  # the same individuals.
  if (n_cate_units > 0L) {
    cate_mean <- colMeans(fits[[1L]]$cate)
    probs <- seq(0, 1, length.out = n_cate_units)
    unit_index <- unique(order(cate_mean)[
      pmax(1L, round(probs * length(cate_mean)))
    ])
    for (chain in seq_along(fits)) {
      for (k in seq_along(unit_index)) {
        label <- sprintf("CATE[unit %d]", unit_index[k])
        monitored[[chain]][[label]] <- fits[[chain]]$cate[, unit_index[k]]
      }
    }
  }

  quantities <- names(monitored[[1L]])
  rows <- lapply(quantities, function(q) {
    chains <- lapply(monitored, function(m) m[[q]])
    zicbcf_diagnose_one(chains, q, geweke_frac1, geweke_frac2)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' Convergence diagnostics for an arbitrary set of monitored quantities
#'
#' The model-agnostic counterpart of \code{\link{zicbcf_convergence}}, for
#' benchmark models whose posterior draws are not returned as a
#' \code{"zicbcf_fit_smear"} object. It reports the same three diagnostics:
#' split-chain Gelman--Rubin \eqn{\widehat{R}}, effective sample size, and the
#' Geweke statistic.
#'
#' @param monitored A named list. Each element is either a numeric vector of
#'   draws from a single chain, or a list of such vectors, one per chain.
#' @param geweke_frac1,geweke_frac2 Fractions of each chain used for the first
#'   and last windows of the Geweke test.
#'
#' @return A data frame in the same format as \code{\link{zicbcf_convergence}}.
#'
#' @export
zicbcf_convergence_table <- function(monitored,
                                     geweke_frac1 = 0.1,
                                     geweke_frac2 = 0.5) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for zicbcf_convergence_table().")
  }
  if (!is.list(monitored) || is.null(names(monitored))) {
    stop("'monitored' must be a named list.")
  }
  rows <- lapply(names(monitored), function(q) {
    chains <- monitored[[q]]
    if (is.numeric(chains)) chains <- list(as.numeric(chains))
    chains <- lapply(chains, as.numeric)
    zicbcf_diagnose_one(chains, q, geweke_frac1, geweke_frac2)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' Scalar functionals monitored by zicbcf_convergence
#'
#' @param fit An object of class \code{"zicbcf_fit_smear"}.
#' @return A named list of numeric vectors, one per retained draw.
#' @keywords internal
zicbcf_monitored_scalars <- function(fit) {
  hurdle_cate <- pnorm(fit$mu_b + fit$tau_b) - pnorm(fit$mu_b)
  list(
    `ATE (response scale)` = fit$ate,
    `Hurdle ATE (probability)` = rowMeans(hurdle_cate),
    `Duan smearing factor` = fit$smearing_factors,
    `sigma_c (log-scale SD)` = as.numeric(fit$sigma_c),
    `Mean mu_b (hurdle prognostic)` = rowMeans(fit$mu_b),
    `Mean tau_b (hurdle treatment)` = rowMeans(fit$tau_b),
    `Mean mu_c (log-scale prognostic)` = rowMeans(fit$mu_c),
    `Mean tau_c (log-scale treatment)` = rowMeans(fit$tau_c)
  )
}


#' Convergence diagnostics for one monitored quantity
#'
#' @param chains List of numeric draw vectors, one per chain.
#' @param quantity Label for the monitored quantity.
#' @param frac1,frac2 Geweke window fractions.
#' @return A one-row data frame.
#' @keywords internal
zicbcf_diagnose_one <- function(chains, quantity, frac1 = 0.1, frac2 = 0.5) {
  n_chains <- length(chains)
  n_draws <- length(chains[[1L]])
  all_draws <- unlist(chains, use.names = FALSE)

  # Split every chain in half so that a single-chain run still yields a
  # split-Rhat, and a multi-chain run yields the split version of the
  # across-chain statistic.
  halves <- unlist(lapply(chains, function(x) {
    mid <- floor(length(x) / 2)
    list(x[seq_len(mid)], x[(mid + 1L):(2L * mid)])
  }), recursive = FALSE)

  constant <- stats::sd(all_draws) < .Machine$double.eps^0.5
  rhat <- NA_real_
  if (!constant) {
    mcmc_halves <- coda::mcmc.list(lapply(halves, coda::mcmc))
    rhat <- tryCatch(
      coda::gelman.diag(mcmc_halves, autoburnin = FALSE,
                        multivariate = FALSE)$psrf[1L, 1L],
      error = function(e) NA_real_
    )
  }

  ess <- if (constant) NA_real_ else {
    sum(vapply(chains, function(x) {
      as.numeric(coda::effectiveSize(coda::mcmc(x)))
    }, numeric(1)))
  }

  geweke <- if (constant) NA_real_ else {
    zs <- vapply(chains, function(x) {
      value <- tryCatch(
        as.numeric(coda::geweke.diag(coda::mcmc(x), frac1 = frac1,
                                     frac2 = frac2)$z),
        error = function(e) NA_real_
      )
      if (length(value) != 1L) NA_real_ else value
    }, numeric(1))
    if (all(is.na(zs))) NA_real_ else max(abs(zs), na.rm = TRUE)
  }

  data.frame(
    quantity = quantity,
    n_chains = n_chains,
    n_draws = n_draws,
    mean = mean(all_draws),
    rhat = rhat,
    ess = ess,
    ess_per_draw = ess / (n_chains * n_draws),
    mcse = if (is.na(ess) || ess <= 0) NA_real_ else stats::sd(all_draws) / sqrt(ess),
    geweke_z = geweke,
    stringsAsFactors = FALSE
  )
}


#' Homoskedasticity diagnostics for Duan's smearing re-transformation
#'
#' Duan's smearing estimator re-transforms log-scale predictions to the response
#' scale with the single scalar \eqn{\widehat{\phi} = n_{+}^{-1}\sum_{i \in
#' \mathcal{A}} \exp(\widehat{\varepsilon}_i)}, one value per posterior draw,
#' applied identically to both counterfactual arms
#' (Duan, 1983). A single scalar is consistent for the
#' response-scale mean when the log-scale errors are independent of the
#' covariates; under heteroskedasticity the correct multiplier varies with
#' \eqn{x}, and every estimate scales linearly in the multiplier that is used.
#'
#' This function tests that assumption and quantifies its consequences. It
#' returns 1. a Breusch--Pagan test in the studentized form of
#' Koenker (1981), which does not itself require
#' normal errors, computed on the posterior-mean log-scale residuals of the
#' active subset, 2. the same test recomputed within every posterior draw, so
#' the evidence is summarized by a posterior distribution of \eqn{p}-values
#' rather than a single plug-in number, and 3. a locally smeared sensitivity
#' analysis that replaces the single scalar with bin-specific smearing factors
#' over strata of the fitted log-scale mean and reports the resulting change in
#' the average treatment effect.
#'
#' The third component is the practically important one. If the locally smeared
#' average treatment effect is close to the scalar-smeared estimate, the
#' homoskedasticity assumption is not materially driving the reported effect,
#' whatever the test says about it.
#'
#' @param fit An object of class \code{"zicbcf_fit_smear"}.
#' @param y The outcome vector supplied to \code{zicbcf_smear()}.
#' @param z The treatment vector supplied to \code{zicbcf_smear()}.
#' @param x Optional design matrix of auxiliary regressors for the
#'   heteroskedasticity test. Defaults to the fitted log-scale mean alone; when
#'   supplied, the test regresses squared residuals on the fitted mean and on
#'   \code{x}, which is the White-style version of the test.
#' @param n_bins Number of strata of the fitted log-scale mean used by the
#'   locally smeared sensitivity analysis.
#' @param n_draws_test Number of posterior draws used for the draw-wise test.
#'   Draws are thinned evenly across the chain.
#'
#' @return A list with elements \code{test} (a one-row data frame with the
#'   posterior-mean-residual Breusch--Pagan statistic, degrees of freedom and
#'   \eqn{p}-value), \code{draw_test} (a data frame summarizing the draw-wise
#'   \eqn{p}-values), \code{residual_scale} (a data frame of residual standard
#'   deviations and smearing factors by stratum of the fitted mean), and
#'   \code{sensitivity} (a data frame comparing the scalar-smeared and locally
#'   smeared average treatment effects).
#'
#' @references
#' Duan, N. (1983). Smearing estimate: a nonparametric retransformation method.
#' \emph{Journal of the American Statistical Association} 78, 605--610.
#'
#' Koenker, R. (1981). A note on studentizing a test for heteroscedasticity.
#' \emph{Journal of Econometrics} 17, 107--112.
#'
#' @export
zicbcf_smearing_diagnostics <- function(fit, y, z, x = NULL,
                                        n_bins = 5L,
                                        n_draws_test = 200L) {
  if (!inherits(fit, "zicbcf_fit_smear")) {
    stop("'fit' must be an object of class 'zicbcf_fit_smear'.")
  }
  y <- as.numeric(y)
  z <- as.numeric(z)
  active <- which(y > 0)
  if (length(active) < 20L) stop("Too few active observations for the diagnostic.")

  log_y <- log(y[active])
  z_active <- z[active]
  n_sim <- nrow(fit$mu_c)

  # Posterior-mean log-scale fitted values and residuals on the active subset.
  fitted_draws <- fit$mu_c[, active, drop = FALSE] +
    matrix(z_active, nrow = n_sim, ncol = length(active), byrow = TRUE) *
    fit$tau_c[, active, drop = FALSE]
  fitted_mean <- colMeans(fitted_draws)
  residual_mean <- log_y - fitted_mean

  test <- zicbcf_breusch_pagan(residual_mean, fitted_mean, x, active)

  # Draw-wise version: one test per retained posterior draw.
  index <- unique(round(seq(1, n_sim, length.out = min(n_draws_test, n_sim))))
  draw_p <- vapply(index, function(s) {
    resid_s <- log_y - fitted_draws[s, ]
    zicbcf_breusch_pagan(resid_s, fitted_draws[s, ], x, active)$p_value
  }, numeric(1))
  draw_test <- data.frame(
    n_draws_tested = length(index),
    median_p = stats::median(draw_p),
    p_below_0.05 = mean(draw_p < 0.05),
    p_below_0.01 = mean(draw_p < 0.01),
    stringsAsFactors = FALSE
  )

  # Residual scale and smearing factor by stratum of the fitted mean.
  breaks <- stats::quantile(fitted_mean, probs = seq(0, 1, length.out = n_bins + 1L))
  bin <- cut(fitted_mean, breaks = unique(breaks), include.lowest = TRUE, labels = FALSE)
  residual_scale <- do.call(rbind, lapply(sort(unique(bin)), function(b) {
    keep <- bin == b
    data.frame(
      stratum = b,
      n = sum(keep),
      fitted_mean_midpoint = stats::median(fitted_mean[keep]),
      residual_sd = stats::sd(residual_mean[keep]),
      smearing_factor = mean(exp(residual_mean[keep])),
      stringsAsFactors = FALSE
    )
  }))

  # Locally smeared sensitivity analysis. Within each posterior draw the scalar
  # factor is replaced by the smearing factor of the stratum each unit's fitted
  # log-scale mean falls into, and the ATE is recomputed on that basis. Values
  # outside the stratum boundaries observed on the active subset are assigned to
  # the nearest terminal stratum, so that counterfactual predictions beyond the
  # observed support still receive the closest available local factor.
  bin_edges <- unique(breaks)
  n_strata <- length(bin_edges) - 1L
  assign_bin <- function(v) {
    clamped <- pmin(pmax(v, bin_edges[1L]), bin_edges[length(bin_edges)])
    as.integer(cut(clamped, breaks = bin_edges, include.lowest = TRUE,
                   labels = FALSE))
  }

  local_ate <- numeric(length(index))
  for (k in seq_along(index)) {
    s <- index[k]
    resid_s <- log_y - fitted_draws[s, ]
    bin_s <- assign_bin(fitted_draws[s, ])

    # One factor per stratum, defaulting to the draw's scalar factor for any
    # stratum that happens to contain no active observation in this draw.
    phi_stratum <- rep(fit$smearing_factors[s], n_strata)
    observed <- tapply(exp(resid_s), factor(bin_s, levels = seq_len(n_strata)), mean)
    phi_stratum[!is.na(observed)] <- observed[!is.na(observed)]

    mu_c_s <- fit$mu_c[s, ]
    tau_c_s <- fit$tau_c[s, ]
    phi0 <- phi_stratum[assign_bin(mu_c_s)]
    phi1 <- phi_stratum[assign_bin(mu_c_s + tau_c_s)]
    p0 <- pnorm(fit$mu_b[s, ])
    p1 <- pnorm(fit$mu_b[s, ] + fit$tau_b[s, ])
    local_ate[k] <- mean(p1 * exp(mu_c_s + tau_c_s) * phi1 -
                           p0 * exp(mu_c_s) * phi0)
  }

  scalar_ate <- fit$ate[index]
  sensitivity <- data.frame(
    scalar_smearing_ate = mean(scalar_ate),
    scalar_smearing_ci_low = unname(stats::quantile(scalar_ate, 0.025)),
    scalar_smearing_ci_high = unname(stats::quantile(scalar_ate, 0.975)),
    local_smearing_ate = mean(local_ate),
    local_smearing_ci_low = unname(stats::quantile(local_ate, 0.025)),
    local_smearing_ci_high = unname(stats::quantile(local_ate, 0.975)),
    relative_change = mean(local_ate) / mean(scalar_ate) - 1,
    stringsAsFactors = FALSE
  )

  list(test = test, draw_test = draw_test, residual_scale = residual_scale,
       sensitivity = sensitivity)
}


#' Studentized Breusch-Pagan test
#'
#' @param residuals Numeric residual vector.
#' @param fitted Numeric fitted-value vector of the same length.
#' @param x Optional matrix of extra regressors (full-sample rows).
#' @param active Row indices of \code{x} corresponding to \code{residuals}.
#' @return A one-row data frame with the statistic, degrees of freedom and
#'   \eqn{p}-value.
#' @keywords internal
zicbcf_breusch_pagan <- function(residuals, fitted, x = NULL, active = NULL) {
  regressors <- cbind(fitted = fitted)
  if (!is.null(x)) {
    x_active <- as.matrix(x)[active, , drop = FALSE]
    keep <- apply(x_active, 2, function(column) stats::sd(column) > 0)
    if (any(keep)) regressors <- cbind(regressors, x_active[, keep, drop = FALSE])
  }
  squared <- residuals^2
  scaled <- squared / mean(squared)
  auxiliary <- stats::lm(scaled ~ regressors)
  r_squared <- summary(auxiliary)$r.squared
  df <- ncol(as.matrix(regressors))
  statistic <- length(residuals) * r_squared
  data.frame(
    statistic = statistic,
    df = df,
    p_value = stats::pchisq(statistic, df = df, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}
