#' Fit Joint Copula-BCF with Selection on Unobservables (Path C)
#'
#' Fits a Heckman-style selection model where the binary selection (hurdle) process
#' and continuous intensity (outcome) process are modeled jointly via a latent bivariate normal distribution
#' using Bayesian Causal Forests.
#'
#' @param y Semicontinuous response variable (must be >= 0)
#' @param z Treatment indicator (binary vector of 0's and 1's)
#' @param x_control Design matrix for the prognostic functions
#' @param x_moderate Design matrix for the treatment effect moderating functions
#' @param pihat_sel Vector of propensity scores for the selection stage (estimated on full sample)
#' @param pihat_out Vector of propensity scores for the outcome stage (estimated on active sample)
#' @param nburn Number of MCMC iterations for burn-in
#' @param nsim Number of MCMC iterations to save after burn-in
#' @param nthin Keep every nthin MCMC draw
#' @param ntree_control Number of trees in prognostic forests
#' @param ntree_moderate Number of trees in treatment moderating forests
#' @param sd_control Prior SD for prognostic forest leaf parameters
#' @param sd_moderate Prior SD for moderating forest leaf parameters
#' @param base_control Base parameter for tree split probability in prognostic forests
#' @param power_control Power parameter for tree split probability in prognostic forests
#' @param base_moderate Base parameter for tree split probability in moderating forests
#' @param power_moderate Power parameter for tree split probability in moderating forests
#' @param nu Degrees of freedom in Inverse-Gamma prior on outcome variance
#' @param lambda Prior scale parameter for outcome variance (estimated from data if NULL)
#' @param sigq Prior quantile for outcome variance calibration
#' @param include_pi "control", "moderate", "both" or "neither" on how to include propensity scores
#' @param update_interval Print status every update_interval iterations
#' @return A list of MCMC draws for Selection and Outcome forests, correlation rho, and variance parameters.
#' @export
pathc_bcf <- function(y, z,
                      x_control, x_moderate = x_control,
                      pihat_sel = NULL, pihat_out = NULL,
                      nburn = 500, nsim = 500, nthin = 1,
                      ntree_control = 250, ntree_moderate = 50,
                      sd_control = 2.0, sd_moderate = 1.0,
                      base_control = 0.95, power_control = 2,
                      base_moderate = 0.25, power_moderate = 3,
                      nu = 3, lambda = NULL, sigq = 0.9,
                      include_pi = "control",
                      update_interval = 100,
                      x_control_est = NULL, x_moderate_est = NULL,
                      z_est = NULL,
                      pihat_sel_est = NULL, pihat_out_est = NULL) {

  n <- length(y)
  if (any(y < 0)) stop("Outcome y must be semicontinuous (>= 0)")
  if (!all(z %in% c(0, 1))) stop("Treatment z must be binary (0 or 1)")

  # Active indicator
  active <- y > 0

  # Default propensity scores if not provided
  if (is.null(pihat_sel)) {
    # Fit simple logistic regression on full sample
    fit_sel <- glm(z ~ x_control, family = binomial(link = "logit"))
    pihat_sel <- predict(fit_sel, type = "response")
  }
  if (is.null(pihat_out)) {
    # For Joint Copula (Path C), the outcome model is fit on the full sample via data augmentation.
    # Therefore, the propensity score used to control for RIC in the outcome prognostic forest
    # must be the full-sample propensity score pihat_sel, not the active-subset propensity score.
    pihat_out <- pihat_sel
  }

  # Estimate out-of-sample propensity scores if prediction is requested
  if (!is.null(x_control_est)) {
    if (is.null(pihat_sel_est)) {
      fit_sel <- glm(z ~ x_control, family = binomial(link = "logit"))
      pihat_sel_est <- predict(fit_sel, newdata = list(x_control = x_control_est), type = "response")
    }
    if (is.null(pihat_out_est)) {
      x_control_active <- x_control[active, , drop = FALSE]
      z_active <- z[active]
      fit_out <- glm(z_active ~ x_control_active, family = binomial(link = "logit"))
      pihat_out_est <- predict(fit_out, newdata = list(x_control_active = x_control_est), type = "response")
    }
  }

  # Build control and moderate design matrices including propensity scores
  x_c_sel <- x_control
  x_m_sel <- x_moderate
  x_c_out <- x_control
  x_m_out <- x_moderate

  if (include_pi %in% c("both", "control")) {
    x_c_sel <- cbind(x_control, pihat_sel)
    x_c_out <- cbind(x_control, pihat_out)
  }
  if (include_pi %in% c("both", "moderate")) {
    x_m_sel <- cbind(x_moderate, pihat_sel)
    x_m_out <- cbind(x_moderate, pihat_out)
  }

  # Build out-of-sample design and basis matrices if provided
  if (!is.null(x_control_est)) {
    if (is.null(x_moderate_est)) x_moderate_est <- x_control_est
    if (is.null(z_est)) stop("z_est (out-of-sample treatment indicator) must be provided when out-of-sample prediction is requested")
    n_est <- nrow(x_control_est)

    x_c_sel_est <- x_control_est
    x_m_sel_est <- x_moderate_est
    x_c_out_est <- x_control_est
    x_m_out_est <- x_moderate_est

    if (include_pi %in% c("both", "control")) {
      x_c_sel_est <- cbind(x_control_est, pihat_sel_est)
      x_c_out_est <- cbind(x_control_est, pihat_out_est)
    }
    if (include_pi %in% c("both", "moderate")) {
      x_m_sel_est <- cbind(x_moderate_est, pihat_sel_est)
      x_m_out_est <- cbind(x_moderate_est, pihat_out_est)
    }

    Omega_sel_con_est <- matrix(rep(1, n_est), nrow = 1)
    Omega_sel_mod_est <- matrix(z_est, nrow = 1)
    Omega_out_con_est <- matrix(rep(1, n_est), nrow = 1)
    Omega_out_mod_est <- matrix(z_est, nrow = 1)

    x_sel_con_est_val <- t(x_c_sel_est)
    x_sel_mod_est_val <- t(x_m_sel_est)
    x_out_con_est_val <- t(x_c_out_est)
    x_out_mod_est_val <- t(x_m_out_est)
  } else {
    n_est <- 0
    Omega_sel_con_est <- matrix(0, nrow = 1, ncol = 0)
    Omega_sel_mod_est <- matrix(0, nrow = 1, ncol = 0)
    Omega_out_con_est <- matrix(0, nrow = 1, ncol = 0)
    Omega_out_mod_est <- matrix(0, nrow = 1, ncol = 0)
    x_sel_con_est_val <- numeric(0)
    x_sel_mod_est_val <- numeric(0)
    x_out_con_est_val <- numeric(0)
    x_out_mod_est_val <- numeric(0)
  }

  # Generate cutpoints lists
  .cp_quantile <- function(x) {
    n <- length(x)
    xs <- sort(unique(x))
    if (length(xs) <= 1) return(xs)
    # Simple quantiles
    q <- quantile(xs, seq(0.01, 0.99, length.out = min(100, length(xs))))
    return(as.numeric(unique(q)))
  }

  cutpoints_sel_con <- lapply(1:ncol(x_c_sel), function(i) .cp_quantile(x_c_sel[, i]))
  cutpoints_sel_mod <- lapply(1:ncol(x_m_sel), function(i) .cp_quantile(x_m_sel[, i]))
  cutpoints_out_con <- lapply(1:ncol(x_c_out), function(i) .cp_quantile(x_c_out[, i]))
  cutpoints_out_mod <- lapply(1:ncol(x_m_out), function(i) .cp_quantile(x_m_out[, i]))

  # We scale the outcome on the log scale for active units
  log_y_active <- log(y[active])
  muy_active <- mean(log_y_active)
  sdy_active <- sd(log_y_active)
  if (is.na(sdy_active) || sdy_active < 1e-6) sdy_active <- 1.0

  # Setup y vector for C++: positive units scaled, zero units preserved
  y_cpp <- y
  y_cpp[active] <- (log_y_active - muy_active) / sdy_active
  y_cpp[!active] <- 0.0

  # Setup basis matrices (Omega)
  # mu forests are vanilla (constant basis = 1)
  # tau forests are modulated (basis = z)
  Omega_sel_con <- matrix(rep(1, n), nrow = 1)
  Omega_sel_mod <- matrix(z, nrow = 1)
  Omega_out_con <- matrix(rep(1, n), nrow = 1)
  Omega_out_mod <- matrix(z, nrow = 1)

  # Prior standard deviation calibration for forests
  # Prognostic prior SD: con_sd^2 / ntree
  # Moderating prior SD: mod_sd^2 / ntree
  Sigma0_sel_con <- matrix(sd_control * sd_control / ntree_control, nrow = 1)
  Sigma0_sel_mod <- matrix(sd_moderate * sd_moderate / ntree_moderate, nrow = 1)
  Sigma0_out_con <- matrix(sd_control * sd_control / ntree_control, nrow = 1)
  Sigma0_out_mod <- matrix(sd_moderate * sd_moderate / ntree_moderate, nrow = 1)

  # Outcome variance (sigma0^2) prior calibration
  if (is.null(lambda)) {
    # Use standard Chipman et al calibration on the positive scaled outcomes
    # Since var(y_cpp[active]) is 1.0, a good prior standard deviation estimate is 1.0
    sighat <- 1.0
    qchi <- qchisq(1.0 - sigq, nu)
    lambda <- (sighat * sighat * qchi) / nu
  }

  # Call C++ core
  res <- pathc_bcfCore(y_ = y_cpp,
                       I_ = as.numeric(active),
                       Omega_sel_con = Omega_sel_con, Omega_sel_mod = Omega_sel_mod,
                       Omega_out_con = Omega_out_con, Omega_out_mod = Omega_out_mod,
                       Omega_sel_con_est = Omega_sel_con_est, Omega_sel_mod_est = Omega_sel_mod_est,
                       Omega_out_con_est = Omega_out_con_est, Omega_out_mod_est = Omega_out_mod_est,
                       x_sel_con_ = t(x_c_sel), x_sel_mod_ = t(x_m_sel),
                       x_out_con_ = t(x_c_out), x_out_mod_ = t(x_m_out),
                       x_sel_con_est_ = x_sel_con_est_val, x_sel_mod_est_ = x_sel_mod_est_val,
                       x_out_con_est_ = x_out_con_est_val, x_out_mod_est_ = x_out_mod_est_val,
                       x_sel_con_info_list = cutpoints_sel_con, x_sel_mod_info_list = cutpoints_sel_mod,
                       x_out_con_info_list = cutpoints_out_con, x_out_mod_info_list = cutpoints_out_mod,
                       burn = nburn, nd = nsim, thin = nthin,
                       ntree_sel_con = ntree_control, ntree_sel_mod = ntree_moderate,
                       ntree_out_con = ntree_control, ntree_out_mod = ntree_moderate,
                       lambda = lambda, nu = nu,
                       Sigma0_sel_con = Sigma0_sel_con, Sigma0_sel_mod = Sigma0_sel_mod,
                       Sigma0_out_con = Sigma0_out_con, Sigma0_out_mod = Sigma0_out_mod,
                       sel_con_alpha = base_control, sel_con_beta = power_control,
                       sel_mod_alpha = base_moderate, sel_mod_beta = power_moderate,
                       out_con_alpha = base_control, out_con_beta = power_control,
                       out_mod_alpha = base_moderate, out_mod_beta = power_moderate,
                       con_scale_df = 1, mod_scale_df = -1,
                       status_interval = update_interval)

  # Re-scale outcome posterior draws back to the original log scale
  res$out_con_post <- muy_active + sdy_active * res$out_con_post
  res$out_mod_post <- sdy_active * res$out_mod_post
  res$yhat_post <- muy_active + sdy_active * res$yhat_post

  # Unmodulated treatment effects
  res$sel_tau_post <- res$sel_tau_post
  res$out_tau_post <- sdy_active * res$out_tau_post

  res$beta_post <- sdy_active * res$beta_post
  res$sigma0_post <- sdy_active * res$sigma0_post
  res$sigma_post <- sdy_active * res$sigma_post

  # For latent augmented outcome V_post, we also re-scale
  res$V_post <- muy_active + sdy_active * res$V_post

  # Propensity scores used
  res$pihat_sel <- pihat_sel
  res$pihat_out <- pihat_out
  res$muy_active <- muy_active
  res$sdy_active <- sdy_active

  if (!is.null(x_control_est)) {
    # Re-scale outcome estimation posterior draws back to the original log scale
    res$out_con_est_post <- muy_active + sdy_active * res$out_con_est_post
    res$out_mod_est_post <- sdy_active * res$out_mod_est_post
    res$yhat_est_post <- res$out_con_est_post + res$out_mod_est_post
    
    res$pihat_sel_est <- pihat_sel_est
    res$pihat_out_est <- pihat_out_est
  }

  return(res)
}
