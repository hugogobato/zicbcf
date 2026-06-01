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

  # 1. Fit the standard ZIC-BCF model using Path A's robust C++ engine
  fit_pathA <- zicbcf_pathA(
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
  nsim_saved <- nrow(fit_pathA$mu_b)
  active_idx <- which(y > 0)
  n_active <- length(active_idx)

  # Hurdle probit predictions
  mu_b <- fit_pathA$mu_b
  tau_b <- fit_pathA$tau_b

  # Continuous log-scale predictions (predicted on the full sample of size n)
  mu_c <- fit_pathA$mu_c
  tau_c <- fit_pathA$tau_c
  sigma_c <- fit_pathA$sigma_c

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
    yhat_hurdle = fit_pathA$yhat_hurdle,
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
    pihat_active = fit_pathA$pihat_active,
    sdy_cont = fit_pathA$sdy_cont,
    muy_cont = fit_pathA$muy_cont,
    raw_fit = fit_pathA$raw_fit
  )

  class(out) <- "zicbcf_fit_smear"
  return(out)
}
