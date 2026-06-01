#' Fit Zero-Inflated Continuous Bayesian Causal Forests (ZIC-BCF)
#'
#' Fits a two-part hurdle Bayesian Causal Forest model where the hurdle
#' is modeled via a probit BCF on the full sample, and the continuous
#' log-intensity is modeled via a Gaussian BCF on the active subset (Y > 0)
#' with Subpopulation Propensity Adjustment (SPA).
#'
#' @param y Semicontinuous response variable (continuous vector >= 0).
#' @param z Binary treatment assignments (0 or 1).
#' @param x_control Design matrix for the prognostic function mu(x).
#' @param x_moderate Design matrix for the treatment effect function tau(x). Defaults to x_control.
#' @param pihat Full-sample propensity score vector. Defaults to 0.5.
#' @param pihat_active Subpopulation propensity score vector on active subset.
#'   If NULL, estimated using logistic regression on the active subset.
#' @param nburn Number of burn-in MCMC iterations.
#' @param nsim Number of MCMC iterations to save after burn-in.
#' @param nthin Thinning interval.
#' @param update_interval Interval to print MCMC status.
#' @param ntree_control_hurdle Number of trees in prognostic hurdle forest.
#' @param ntree_moderate_hurdle Number of trees in moderating hurdle forest.
#' @param ntree_control_continuous Number of trees in prognostic continuous forest.
#' @param ntree_moderate_continuous Number of trees in moderating continuous forest.
#' @param sd_control_hurdle Prior median of SD(mu_b(x)) on the latent probit scale.
#' @param sd_moderate_hurdle Prior median of SD(tau_b(x)) on the latent probit scale.
#' @param sd_control_continuous Prior median of SD(mu_c(x)) on the log-outcome scale.
#'   Defaults to 2 * sd(log(y[y > 0])).
#' @param sd_moderate_continuous Prior median of SD(tau_c(x)) on the log-outcome scale.
#'   Defaults to 0.25 * sd(log(y[y > 0])) / sd(z[y > 0]).
#' @param con_alpha_hurdle Base for tree prior on prognostic hurdle forest.
#' @param con_beta_hurdle Power for tree prior on prognostic hurdle forest.
#' @param mod_alpha_hurdle Base for tree prior on moderating hurdle forest.
#' @param mod_beta_hurdle Power for tree prior on moderating hurdle forest.
#' @param con_alpha_continuous Base for tree prior on prognostic continuous forest.
#' @param con_beta_continuous Power for tree prior on prognostic continuous forest.
#' @param mod_alpha_continuous Base for tree prior on moderating continuous forest.
#' @param mod_beta_continuous Power for tree prior on moderating continuous forest.
#' @param nu Degrees of freedom in continuous variance prior.
#' @param lambda Scale parameter in continuous variance prior. Calibrated if NULL.
#' @param sigq Calibration quantile for continuous variance prior.
#' @param sighat Calibration estimate for continuous variance.
#' @param include_pi "control", "moderate", "both" or "none".
#' @param use_muscale Use half-Cauchy hyperprior on the scale of mu.
#' @param use_tauscale Use half-Normal hyperprior on the scale of tau.
#'
#' @return A list of class \code{"zicbcf_fit"} with draws.
#' @export
zicbcf_fit <- function(
    y, z,
    x_control, x_moderate = x_control,
    pihat = rep(0.5, length(y)),
    pihat_active = NULL,
    nburn = 500, nsim = 1000, nthin = 1,
    update_interval = 100,
    ntree_control_hurdle = 200, ntree_moderate_hurdle = 50,
    ntree_control_continuous = 200, ntree_moderate_continuous = 50,
    sd_control_hurdle = 2.0, sd_moderate_hurdle = 1.0,
    sd_control_continuous = NULL, sd_moderate_continuous = NULL,
    con_alpha_hurdle = 0.95, con_beta_hurdle = 2.0,
    mod_alpha_hurdle = 0.25, mod_beta_hurdle = 3.0,
    con_alpha_continuous = 0.95, con_beta_continuous = 2.0,
    mod_alpha_continuous = 0.25, mod_beta_continuous = 3.0,
    nu = 3, lambda = NULL, sigq = 0.9, sighat = NULL,
    include_pi = "control",
    use_muscale = TRUE, use_tauscale = TRUE
) {

  # -- input validation --------------------------------------------------------
  x_control <- as.matrix(x_control)
  x_moderate <- as.matrix(x_moderate)
  pihat <- as.matrix(pihat)

  n <- length(y)
  if (!.ident(n, nrow(x_control), nrow(x_moderate), nrow(pihat))) {
    stop("Data size mismatch. length(y), nrow(x_control), nrow(x_moderate), nrow(pihat) must all match.")
  }
  if (length(z) != n) stop("length(z) must equal length(y)")

  if (any(is.na(y))) stop("Missing values in y")
  if (any(is.na(z))) stop("Missing values in z")
  if (any(is.na(x_control))) stop("Missing values in x_control")
  if (any(is.na(x_moderate))) stop("Missing values in x_moderate")
  if (any(is.na(pihat))) stop("Missing values in pihat")

  if (any(!is.finite(y))) stop("Non-finite values in y")
  if (any(!is.finite(z))) stop("Non-finite values in z")
  if (any(!is.finite(x_control))) stop("Non-finite values in x_control")
  if (any(!is.finite(x_moderate))) stop("Non-finite values in x_moderate")
  if (any(!is.finite(pihat))) stop("Non-finite values in pihat")

  if (any(y < 0)) stop("Outcome y must be non-negative")
  if (!all(z %in% c(0, 1))) stop("Treatment z must contain only 0 and 1")
  if (nburn < 0) stop("nburn must be non-negative")
  if (nsim < 1) stop("nsim must be positive")
  if (nthin < 1) stop("nthin must be positive")
  if (nburn < 100) warning("A low (<100) value for nburn was supplied")

  # -- active subset & SPA propensity -----------------------------------------
  active_idx <- which(y > 0)
  n_c <- length(active_idx)
  if (n_c < 10) stop("Too few positive outcomes (<10) to fit continuous intensity stage.")

  x_df <- as.data.frame(x_control)
  names(x_df) <- paste0("V", seq_len(ncol(x_df)))

  if (is.null(pihat_active)) {
    x_active_df <- x_df[active_idx, , drop = FALSE]
    z_active <- z[active_idx]

    fit_spa <- tryCatch(
      suppressWarnings(glm(z_active ~ ., data = x_active_df, family = binomial(link = "logit"))),
      error = function(e) NULL
    )

    if (!is.null(fit_spa)) {
      pihat_active <- suppressWarnings(predict(fit_spa, type = "response"))
      pihat_active[is.na(pihat_active)] <- mean(z_active)
      pihat_active <- pmin(pmax(pihat_active, 0.01), 0.99)
    } else {
      pihat_active <- rep(mean(z_active), n_c)
    }
  } else {
    if (length(pihat_active) == n) {
      pihat_active <- pihat_active[active_idx]
    }
    if (length(pihat_active) != n_c) {
      stop("pihat_active must have length equal to the active subset or full sample")
    }
  }
  pihat_active <- as.matrix(pihat_active)

  # Active-subset propensity model predicted for ALL observations (for out-of-sample est)
  fit_spa_full <- tryCatch(
    suppressWarnings(glm(z ~ ., data = x_df, subset = active_idx, family = binomial(link = "logit"))),
    error = function(e) NULL
  )

  if (!is.null(fit_spa_full)) {
    pihat_active_full <- suppressWarnings(predict(fit_spa_full, newdata = x_df, type = "response"))
    pihat_active_full[is.na(pihat_active_full)] <- mean(z[active_idx])
    pihat_active_full <- pmin(pmax(pihat_active_full, 0.01), 0.99)
  } else {
    pihat_active_full <- rep(mean(z[active_idx]), n)
  }
  pihat_active_full <- as.matrix(pihat_active_full)

  # Disable pi inclusion if degenerate
  if (include_pi != "none" && length(unique(pihat)) == 1) {
    warning("All values of pihat are equal. pihat will not be included among covariates.")
    include_pi <- "none"
  }

  # -- Hurdle designs ----------------------------------------------------------
  x_c_hurdle <- x_control
  x_m_hurdle <- x_moderate
  if (include_pi %in% c("both", "control")) x_c_hurdle <- cbind(x_control, pihat)
  if (include_pi %in% c("both", "moderate")) x_m_hurdle <- cbind(x_moderate, pihat)

  cutpoint_list_c_hurdle <- lapply(seq_len(ncol(x_c_hurdle)), function(i) .cp_quantile(x_c_hurdle[, i]))
  cutpoint_list_m_hurdle <- lapply(seq_len(ncol(x_m_hurdle)), function(i) .cp_quantile(x_m_hurdle[, i]))

  Omega_con_hurdle <- matrix(rep(1, n), ncol = 1)
  Omega_mod_hurdle <- matrix(z, ncol = 1)

  # Half-Normal prior median calibration:
  # For HN(sigma), median = 0.6745 * sigma. To make sd_moderate the prior MEDIAN of
  # |tau scale|, set sigma = sd_moderate / 0.6745 when use_tauscale = TRUE.
  hn_factor <- if (use_tauscale) 0.6745 else 1.0
  sd_mod_hurdle_eff <- sd_moderate_hurdle / hn_factor
  Sigma0_con_hurdle <- matrix(sd_control_hurdle^2 / ntree_control_hurdle, nrow = 1)
  Sigma0_mod_hurdle <- matrix(sd_mod_hurdle_eff^2 / ntree_moderate_hurdle, nrow = 1)

  # -- Continuous designs (active subset) --------------------------------------
  y_cont <- log(y[active_idx])
  if (!all(is.finite(y_cont))) stop("Non-finite log(y) in active subset.")
  muy_cont <- mean(y_cont)
  sdy_cont <- sd(y_cont)
  if (is.na(sdy_cont) || sdy_cont < 1e-8) sdy_cont <- 1.0
  y_cont_scale <- (y_cont - muy_cont) / sdy_cont

  # Data-adaptive defaults for continuous-stage scales (computed on log-scale)
  if (is.null(sd_control_continuous)) sd_control_continuous <- 2.0 * sdy_cont
  if (is.null(sd_moderate_continuous)) {
    z_act <- z[active_idx]
    sdz_act <- sd(z_act)
    if (is.na(sdz_act) || sdz_act < 1e-8) sdz_act <- 0.5
    sd_moderate_continuous <- 0.25 * sdy_cont / sdz_act
  }

  x_c_cont <- x_control[active_idx, , drop = FALSE]
  x_m_cont <- x_moderate[active_idx, , drop = FALSE]
  if (include_pi %in% c("both", "control")) x_c_cont <- cbind(x_c_cont, pihat_active)
  if (include_pi %in% c("both", "moderate")) x_m_cont <- cbind(x_m_cont, pihat_active)

  cutpoint_list_c_cont <- lapply(seq_len(ncol(x_c_cont)), function(i) .cp_quantile(x_c_cont[, i]))
  cutpoint_list_m_cont <- lapply(seq_len(ncol(x_m_cont)), function(i) .cp_quantile(x_m_cont[, i]))

  Omega_con_cont <- matrix(rep(1, n_c), ncol = 1)
  Omega_mod_cont <- matrix(z[active_idx], ncol = 1)

  con_sd_cont <- sd_control_continuous / sdy_cont
  mod_sd_cont <- (sd_moderate_continuous / sdy_cont) / hn_factor

  Sigma0_con_cont <- matrix(con_sd_cont^2 / ntree_control_continuous, nrow = 1)
  Sigma0_mod_cont <- matrix(mod_sd_cont^2 / ntree_moderate_continuous, nrow = 1)

  # -- Continuous out-of-sample design (all n observations) --------------------
  x_c_cont_est <- x_control
  x_m_cont_est <- x_moderate
  if (include_pi %in% c("both", "control")) x_c_cont_est <- cbind(x_control, pihat_active_full)
  if (include_pi %in% c("both", "moderate")) x_m_cont_est <- cbind(x_moderate, pihat_active_full)

  Omega_con_cont_est <- matrix(rep(1, n), ncol = 1)
  Omega_mod_cont_est <- matrix(z, ncol = 1)

  # -- Sigma prior calibration (Chipman-style) ---------------------------------
  if (is.null(lambda)) {
    if (is.null(sighat)) {
      lmf <- tryCatch(
        suppressWarnings(lm(y_cont_scale ~ x_c_cont)),
        error = function(e) NULL
      )
      sighat_lm <- if (!is.null(lmf)) {
        s <- summary(lmf)$sigma
        if (is.na(s) || s <= 0) 0.9 else s
      } else {
        0.9
      }
      sighat <- sighat_lm * sdy_cont
    }
    qchi <- qchisq(1.0 - sigq, nu)
    lambda <- ((sighat / sdy_cont)^2 * qchi) / nu
  }

  y_hurdle_indicator <- as.numeric(y > 0)

  # -- Call C++ core -----------------------------------------------------------
  fit <- zicbcfCore(
    y_hurdle = y_hurdle_indicator,
    Omega_con_hurdle = t(Omega_con_hurdle),
    Omega_mod_hurdle = t(Omega_mod_hurdle),
    x_con_hurdle_ = t(x_c_hurdle),
    x_mod_hurdle_ = t(x_m_hurdle),
    x_con_info_hurdle_list = cutpoint_list_c_hurdle,
    x_mod_info_hurdle_list = cutpoint_list_m_hurdle,
    ntree_con_hurdle = ntree_control_hurdle,
    ntree_mod_hurdle = ntree_moderate_hurdle,
    Sigma0_con_hurdle = Sigma0_con_hurdle,
    Sigma0_mod_hurdle = Sigma0_mod_hurdle,
    con_alpha_hurdle = con_alpha_hurdle,
    con_beta_hurdle = con_beta_hurdle,
    mod_alpha_hurdle = mod_alpha_hurdle,
    mod_beta_hurdle = mod_beta_hurdle,
    vanilla_hurdle = TRUE,
    use_con_scale_hurdle = use_muscale,
    use_mod_scale_hurdle = use_tauscale,
    con_scale_df_hurdle = 1,
    mod_scale_df_hurdle = -1,

    y_continuous = y_cont_scale,
    Omega_con_continuous = t(Omega_con_cont),
    Omega_mod_continuous = t(Omega_mod_cont),
    x_con_continuous_ = t(x_c_cont),
    x_mod_continuous_ = t(x_m_cont),
    x_con_info_continuous_list = cutpoint_list_c_cont,
    x_mod_info_continuous_list = cutpoint_list_m_cont,
    ntree_con_continuous = ntree_control_continuous,
    ntree_mod_continuous = ntree_moderate_continuous,
    Sigma0_con_continuous = Sigma0_con_cont,
    Sigma0_mod_continuous = Sigma0_mod_cont,
    con_alpha_continuous = con_alpha_continuous,
    con_beta_continuous = con_beta_continuous,
    mod_alpha_continuous = mod_alpha_continuous,
    mod_beta_continuous = mod_beta_continuous,
    vanilla_continuous = TRUE,
    use_con_scale_continuous = use_muscale,
    use_mod_scale_continuous = use_tauscale,
    con_scale_df_continuous = 1,
    mod_scale_df_continuous = -1,

    x_con_continuous_est_ = t(x_c_cont_est),
    x_mod_continuous_est_ = t(x_m_cont_est),
    Omega_con_continuous_est = t(Omega_con_cont_est),
    Omega_mod_continuous_est = t(Omega_mod_cont_est),

    burn = nburn,
    nd = nsim,
    thin = nthin,
    lambda = lambda,
    nu = nu,
    status_interval = update_interval
  )

  # -- Post-processing: reconstruct potential outcomes -------------------------
  mu_b <- fit$m_hurdle_post
  tau_b <- fit$b_hurdle_post

  mu_c <- muy_cont + sdy_cont * fit$m_continuous_est_post
  tau_c <- sdy_cont * fit$b_continuous_est_post
  sigma_c <- sdy_cont * fit$sigma_continuous

  p0 <- pnorm(mu_b)
  p1 <- pnorm(mu_b + tau_b)
  half_var <- 0.5 * sigma_c^2
  y0_plus <- exp(mu_c + half_var)
  y1_plus <- exp(mu_c + tau_c + half_var)

  mu0 <- p0 * y0_plus
  mu1 <- p1 * y1_plus
  cate <- mu1 - mu0
  ate <- rowMeans(cate)

  out <- list(
    yhat_hurdle = fit$yhat_hurdle_post,
    mu_b = mu_b,
    tau_b = tau_b,
    mu_c = mu_c,
    tau_c = tau_c,
    sigma_c = sigma_c,
    mu0 = mu0,
    mu1 = mu1,
    cate = cate,
    ate = ate,
    pihat_active = pihat_active_full,
    sdy_cont = sdy_cont,
    muy_cont = muy_cont,
    raw_fit = fit
  )

  class(out) <- "zicbcf_fit"
  return(out)
}

#' Fit Zero-Inflated Continuous BCF with Duan's Smearing (ZIC-BCF-Smear)
#'
#' Fits a two-part hurdle Bayesian Causal Forest model where the hurdle
#' is modeled via a probit BCF on the full sample, and the continuous
#' log-intensity is modeled via a Gaussian BCF on the active subset (Y > 0)
#' with Subpopulation Propensity Adjustment (SPA). The outcome is re-transformed
#' to the original scale using Duan's non-parametric smearing estimator to avoid
#' distributional misspecification bias (e.g. under Tweedie or Gamma).
#'
#' @param y Semicontinuous response variable (continuous vector >= 0).
#' @param z Binary treatment assignments (0 or 1).
#' @param x_control Design matrix for the prognostic function mu(x).
#' @param x_moderate Design matrix for the treatment effect function tau(x). Defaults to x_control.
#' @param pihat Full-sample propensity score vector. Defaults to 0.5.
#' @param pihat_active Subpopulation propensity score vector on active subset.
#'   If NULL, estimated using logistic regression on the active subset.
#' @param nburn Number of burn-in MCMC iterations.
#' @param nsim Number of MCMC iterations to save after burn-in.
#' @param nthin Thinning interval.
#' @param update_interval Interval to print MCMC status.
#' @param ntree_control_hurdle Number of trees in prognostic hurdle forest.
#' @param ntree_moderate_hurdle Number of trees in moderating hurdle forest.
#' @param ntree_control_continuous Number of trees in prognostic continuous forest.
#' @param ntree_moderate_continuous Number of trees in moderating continuous forest.
#' @param sd_control_hurdle Prior median of SD(mu_b(x)) on the latent probit scale.
#' @param sd_moderate_hurdle Prior median of SD(tau_b(x)) on the latent probit scale.
#' @param sd_control_continuous Prior median of SD(mu_c(x)) on the log-outcome scale.
#'   Defaults to 2 * sd(log(y[y > 0])).
#' @param sd_moderate_continuous Prior median of SD(tau_c(x)) on the log-outcome scale.
#'   Defaults to 0.25 * sd(log(y[y > 0])) / sd(z[y > 0]).
#' @param con_alpha_hurdle Base for tree prior on prognostic hurdle forest.
#' @param con_beta_hurdle Power for tree prior on prognostic hurdle forest.
#' @param mod_alpha_hurdle Base for tree prior on moderating hurdle forest.
#' @param mod_beta_hurdle Power for tree prior on moderating hurdle forest.
#' @param con_alpha_continuous Base for tree prior on prognostic continuous forest.
#' @param con_beta_continuous Power for tree prior on prognostic continuous forest.
#' @param mod_alpha_continuous Base for tree prior on moderating continuous forest.
#' @param mod_beta_continuous Power for tree prior on moderating continuous forest.
#' @param nu Degrees of freedom in continuous variance prior.
#' @param lambda Scale parameter in continuous variance prior. Calibrated if NULL.
#' @param sigq Calibration quantile for continuous variance prior.
#' @param sighat Calibration estimate for continuous variance.
#' @param include_pi "control", "moderate", "both" or "none".
#' @param use_muscale Use half-Cauchy hyperprior on the scale of mu.
#' @param use_tauscale Use half-Normal hyperprior on the scale of tau.
#'
#' @return A list of class \code{"zicbcf_fit_smear"} with draws.
#' @export
zicbcf_smear <- function(
    y, z,
    x_control, x_moderate = x_control,
    pihat = rep(0.5, length(y)),
    pihat_active = NULL,
    nburn = 500, nsim = 1000, nthin = 1,
    update_interval = 100,
    ntree_control_hurdle = 200, ntree_moderate_hurdle = 50,
    ntree_control_continuous = 200, ntree_moderate_continuous = 50,
    sd_control_hurdle = 2.0, sd_moderate_hurdle = 1.0,
    sd_control_continuous = NULL, sd_moderate_continuous = NULL,
    con_alpha_hurdle = 0.95, con_beta_hurdle = 2.0,
    mod_alpha_hurdle = 0.25, mod_beta_hurdle = 3.0,
    con_alpha_continuous = 0.95, con_beta_continuous = 2.0,
    mod_alpha_continuous = 0.25, mod_beta_continuous = 3.0,
    nu = 3, lambda = NULL, sigq = 0.9, sighat = NULL,
    include_pi = "control",
    use_muscale = TRUE, use_tauscale = TRUE
) {

  # 1. Fit the standard ZIC-BCF model using our robust C++ engine
  fit_raw <- zicbcf_fit(
    y = y, z = z,
    x_control = x_control, x_moderate = x_moderate,
    pihat = pihat, pihat_active = pihat_active,
    nburn = nburn, nsim = nsim, nthin = nthin,
    update_interval = update_interval,
    ntree_control_hurdle = ntree_control_hurdle,
    ntree_moderate_hurdle = ntree_moderate_hurdle,
    ntree_control_continuous = ntree_control_continuous,
    ntree_moderate_continuous = ntree_moderate_continuous,
    sd_control_hurdle = sd_control_hurdle,
    sd_moderate_hurdle = sd_moderate_hurdle,
    sd_control_continuous = sd_control_continuous,
    sd_moderate_continuous = sd_moderate_continuous,
    con_alpha_hurdle = con_alpha_hurdle,
    con_beta_hurdle = con_beta_hurdle,
    mod_alpha_hurdle = mod_alpha_hurdle,
    mod_beta_hurdle = mod_beta_hurdle,
    con_alpha_continuous = con_alpha_continuous,
    con_beta_continuous = con_beta_continuous,
    mod_alpha_continuous = mod_alpha_continuous,
    mod_beta_continuous = mod_beta_continuous,
    nu = nu, lambda = lambda, sigq = sigq, sighat = sighat,
    include_pi = include_pi,
    use_muscale = use_muscale, use_tauscale = use_tauscale
  )

  # 2. Extract indices and sizes
  n <- length(y)
  nsim_saved <- nrow(fit_raw$mu_b)
  active_idx <- which(y > 0)
  n_active <- length(active_idx)

  # Hurdle probit predictions
  mu_b <- fit_raw$mu_b
  tau_b <- fit_raw$tau_b

  # Continuous log-scale predictions (predicted on the full sample of size n)
  mu_c <- fit_raw$mu_c
  tau_c <- fit_raw$tau_c
  sigma_c <- fit_raw$sigma_c

  # 3. Retrieve observed active outcomes
  log_y_active <- log(y[active_idx])

  # 4. Perform Duan's Non-Parametric Smearing Re-transformation
  mu0 <- matrix(0, nrow = nsim_saved, ncol = n)
  mu1 <- matrix(0, nrow = nsim_saved, ncol = n)
  cate <- matrix(0, nrow = nsim_saved, ncol = n)
  smearing_factors <- rep(0, nsim_saved)

  for (s in 1:nsim_saved) {
    # Predicted log-outcome on the active subset under actual treatment z
    z_active <- z[active_idx]
    pred_log_active <- mu_c[s, active_idx] + z_active * tau_c[s, active_idx]

    # Compute active log-scale residuals
    residuals_active <- log_y_active - pred_log_active

    # Compute Duan's smearing factor for this iteration
    phi_smear <- mean(exp(residuals_active))
    smearing_factors[s] <- phi_smear

    # Re-transform potential outcomes to response scale using non-parametric smearing
    p0 <- pnorm(mu_b[s, ])
    p1 <- pnorm(mu_b[s, ] + tau_b[s, ])

    y0_plus <- exp(mu_c[s, ]) * phi_smear
    y1_plus <- exp(mu_c[s, ] + tau_c[s, ]) * phi_smear

    mu0[s, ] <- p0 * y0_plus
    mu1[s, ] <- p1 * y1_plus
    cate[s, ] <- mu1[s, ] - mu0[s, ]
  }

  ate <- rowMeans(cate)

  out <- list(
    yhat_hurdle = fit_raw$yhat_hurdle,
    mu_b = mu_b,
    tau_b = tau_b,
    mu_c = mu_c,
    tau_c = tau_c,
    sigma_c = sigma_c,
    smearing_factors = smearing_factors,
    mu0 = mu0,
    mu1 = mu1,
    cate = cate,
    ate = ate,
    pihat_active = fit_raw$pihat_active,
    sdy_cont = fit_raw$sdy_cont,
    muy_cont = fit_raw$muy_cont,
    raw_fit = fit_raw$raw_fit
  )

  class(out) <- "zicbcf_fit_smear"
  return(out)
}
