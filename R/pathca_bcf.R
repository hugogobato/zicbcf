#' Fit Joint Copula-BCF combining Path E's active-outcome collapse with
#' Path A's variance-control priors (Path C-A)
#'
#' Path C-A extends Path E (Collapsed Joint Copula-BCF with active-subset
#' outcome forests) by adopting three variance-reduction ingredients from
#' Path A:
#' \enumerate{
#'   \item The conjugate update for the joint covariance (beta, sigma0^2)
#'         is restricted to the active subset (Y > 0). Augmented V's for
#'         inactive units carry no independent information about beta and
#'         only reinforce its current value, amplifying noise under model
#'         misspecification (e.g., Tweedie).
#'   \item Data-adaptive prior calibration on the log-outcome scale,
#'         matching Path A: \code{sd_control_continuous = 2 * sd(log(Y[Y>0]))},
#'         \code{sd_moderate_continuous = 0.25 * sd(log(Y[Y>0])) / sd(Z[Y>0])},
#'         and Chipman-style \code{lambda} calibration from OLS on the active
#'         subset.
#'   \item Subpopulation Propensity Adjustment (SPA): the outcome forests use
#'         a propensity score estimated on the active subset
#'         (\code{pihat_active}), while the selection forests use the
#'         full-sample propensity (\code{pihat_sel}).
#' }
#'
#' @param y Semicontinuous response variable (must be >= 0)
#' @param z Treatment indicator (binary vector of 0's and 1's)
#' @param x_control Design matrix for the prognostic functions
#' @param x_moderate Design matrix for the treatment effect moderating functions
#' @param pihat_sel Full-sample propensity score (used by selection forests).
#'   Estimated by logistic regression on the full sample if NULL.
#' @param pihat_active Active-subset propensity score (used by outcome forests
#'   under SPA). Estimated on the active subset if NULL. May be a length-n
#'   vector or length-n_active vector.
#' @param nburn Number of MCMC iterations for burn-in
#' @param nsim Number of MCMC iterations to save after burn-in
#' @param nthin Keep every nthin MCMC draw
#' @param ntree_control Number of trees in prognostic forests
#' @param ntree_moderate Number of trees in treatment moderating forests
#' @param sd_control_hurdle Prior SD for the selection prognostic forest
#' @param sd_moderate_hurdle Prior SD for the selection moderating forest
#' @param sd_control_continuous Prior SD for the outcome prognostic forest on
#'   the log-outcome scale. Defaults to \code{2 * sd(log(y[y > 0]))}.
#' @param sd_moderate_continuous Prior SD for the outcome moderating forest
#'   on the log-outcome scale. Defaults to
#'   \code{0.25 * sd(log(y[y > 0])) / sd(z[y > 0])}.
#' @param base_control Base for tree-prior on prognostic forests
#' @param power_control Power for tree-prior on prognostic forests
#' @param base_moderate Base for tree-prior on moderating forests
#' @param power_moderate Power for tree-prior on moderating forests
#' @param nu Degrees of freedom in Inverse-Gamma prior on outcome variance
#' @param lambda Prior scale for outcome variance (Chipman-calibrated from
#'   OLS on the active subset if NULL)
#' @param sigq Prior quantile for outcome variance calibration
#' @param include_pi One of "control", "moderate", "both", "neither" — how to
#'   include propensity scores in the design matrices.
#' @param update_interval Print status every update_interval iterations
#' @return A list of MCMC draws for Selection and Outcome forests, correlation
#'   rho, and variance parameters.
#' @export
pathca_bcf <- function(y, z,
                       x_control, x_moderate = x_control,
                       pihat_sel = NULL,
                       pihat_active = NULL,
                       nburn = 500, nsim = 500, nthin = 1,
                       ntree_control = 250, ntree_moderate = 50,
                       sd_control_hurdle = 2.0, sd_moderate_hurdle = 1.0,
                       sd_control_continuous = NULL,
                       sd_moderate_continuous = NULL,
                       base_control = 0.95, power_control = 2,
                       base_moderate = 0.25, power_moderate = 3,
                       nu = 3, lambda = NULL, sigq = 0.9,
                       include_pi = "control",
                       update_interval = 100,
                       x_control_est = NULL, x_moderate_est = NULL,
                       z_est = NULL,
                       pihat_sel_est = NULL, pihat_active_est = NULL) {

  n <- length(y)
  if (any(y < 0)) stop("Outcome y must be semicontinuous (>= 0)")
  if (!all(z %in% c(0, 1))) stop("Treatment z must be binary (0 or 1)")

  active <- y > 0
  n_active <- sum(active)
  if (n_active < 10) stop("Too few positive outcomes (<10) to fit active outcome stage.")

  # --- Full-sample propensity (for selection forests) ------------------------
  if (is.null(pihat_sel)) {
    fit_sel <- suppressWarnings(
      glm(z ~ x_control, family = binomial(link = "logit"))
    )
    pihat_sel <- predict(fit_sel, type = "response")
  }

  # --- SPA (active-subset) propensity (for outcome forests) ------------------
  x_df <- as.data.frame(x_control)
  names(x_df) <- paste0("V", seq_len(ncol(x_df)))

  if (is.null(pihat_active)) {
    x_active_df <- x_df[active, , drop = FALSE]
    z_active <- z[active]
    fit_spa_full <- tryCatch(
      suppressWarnings(glm(z ~ ., data = x_df, subset = which(active),
                           family = binomial(link = "logit"))),
      error = function(e) NULL
    )
    if (!is.null(fit_spa_full)) {
      pihat_active <- suppressWarnings(predict(fit_spa_full, newdata = x_df, type = "response"))
      pihat_active[is.na(pihat_active)] <- mean(z_active)
      pihat_active <- pmin(pmax(pihat_active, 0.01), 0.99)
    } else {
      pihat_active <- rep(mean(z_active), n)
    }
  } else if (length(pihat_active) == n_active) {
    p_full <- numeric(n)
    p_full[active] <- pihat_active
    p_full[!active] <- mean(z[active])
    pihat_active <- p_full
  } else if (length(pihat_active) != n) {
    stop("pihat_active must have length n or n_active")
  }

  # Disable pi if degenerate
  if (include_pi != "neither" && length(unique(pihat_sel)) == 1) {
    warning("All values of pihat_sel are equal. pihat will not be included among covariates.")
    include_pi <- "neither"
  }

  # --- Design matrices: selection uses pihat_sel, outcome uses pihat_active --
  x_c_sel <- x_control
  x_m_sel <- x_moderate
  x_c_out_full <- x_control
  x_m_out_full <- x_moderate

  if (include_pi %in% c("both", "control")) {
    x_c_sel <- cbind(x_control, pihat_sel)
    x_c_out_full <- cbind(x_control, pihat_active)
  }
  if (include_pi %in% c("both", "moderate")) {
    x_m_sel <- cbind(x_moderate, pihat_sel)
    x_m_out_full <- cbind(x_moderate, pihat_active)
  }

  x_c_out_active <- x_c_out_full[active, , drop = FALSE]
  x_m_out_active <- x_m_out_full[active, , drop = FALSE]

  # --- Out-of-sample design ---------------------------------------------------
  if (!is.null(x_control_est)) {
    if (is.null(x_moderate_est)) x_moderate_est <- x_control_est
    if (is.null(z_est)) stop("z_est required for out-of-sample prediction")
    n_est <- nrow(x_control_est)

    if (is.null(pihat_sel_est)) {
      fit_sel_est <- suppressWarnings(
        glm(z ~ x_control, family = binomial(link = "logit"))
      )
      pihat_sel_est <- predict(fit_sel_est,
                               newdata = list(x_control = x_control_est),
                               type = "response")
    }

    if (is.null(pihat_active_est)) {
      x_active_df <- x_df[active, , drop = FALSE]
      z_active <- z[active]
      fit_spa_est <- tryCatch(
        suppressWarnings(glm(z_active ~ ., data = x_active_df,
                             family = binomial(link = "logit"))),
        error = function(e) NULL
      )
      if (!is.null(fit_spa_est)) {
        x_est_df <- as.data.frame(x_control_est)
        names(x_est_df) <- names(x_active_df)
        pihat_active_est <- suppressWarnings(
          predict(fit_spa_est, newdata = x_est_df, type = "response")
        )
        pihat_active_est[is.na(pihat_active_est)] <- mean(z_active)
        pihat_active_est <- pmin(pmax(pihat_active_est, 0.01), 0.99)
      } else {
        pihat_active_est <- rep(mean(z[active]), n_est)
      }
    }

    x_c_sel_est <- x_control_est
    x_m_sel_est <- x_moderate_est
    x_c_out_est <- x_control_est
    x_m_out_est <- x_moderate_est

    if (include_pi %in% c("both", "control")) {
      x_c_sel_est <- cbind(x_control_est, pihat_sel_est)
      x_c_out_est <- cbind(x_control_est, pihat_active_est)
    }
    if (include_pi %in% c("both", "moderate")) {
      x_m_sel_est <- cbind(x_moderate_est, pihat_sel_est)
      x_m_out_est <- cbind(x_moderate_est, pihat_active_est)
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

  # --- Cutpoints --------------------------------------------------------------
  .cp_quantile <- function(x) {
    xs <- sort(unique(x))
    if (length(xs) <= 1) return(xs)
    q <- quantile(xs, seq(0.01, 0.99, length.out = min(100, length(xs))))
    return(as.numeric(unique(q)))
  }

  cutpoints_sel_con <- lapply(1:ncol(x_c_sel), function(i) .cp_quantile(x_c_sel[, i]))
  cutpoints_sel_mod <- lapply(1:ncol(x_m_sel), function(i) .cp_quantile(x_m_sel[, i]))
  cutpoints_out_con <- lapply(1:ncol(x_c_out_full), function(i) .cp_quantile(x_c_out_full[, i]))
  cutpoints_out_mod <- lapply(1:ncol(x_m_out_full), function(i) .cp_quantile(x_m_out_full[, i]))

  # --- Scale outcome on log scale (active units) ------------------------------
  log_y_active <- log(y[active])
  muy_active <- mean(log_y_active)
  sdy_active <- sd(log_y_active)
  if (is.na(sdy_active) || sdy_active < 1e-6) sdy_active <- 1.0

  y_cpp <- y
  y_cpp[active] <- (log_y_active - muy_active) / sdy_active
  y_cpp[!active] <- 0.0

  # --- Path A-style data-adaptive prior scales --------------------------------
  if (is.null(sd_control_continuous)) sd_control_continuous <- 2.0 * sdy_active
  if (is.null(sd_moderate_continuous)) {
    z_act <- z[active]
    sdz_act <- sd(z_act)
    if (is.na(sdz_act) || sdz_act < 1e-8) sdz_act <- 0.5
    sd_moderate_continuous <- 0.25 * sdy_active / sdz_act
  }

  # Convert outcome prior SDs to scaled-outcome units (the C++ core operates
  # on (log y - muy)/sdy with target std-dev close to 1).
  con_sd_out <- sd_control_continuous / sdy_active
  mod_sd_out <- sd_moderate_continuous / sdy_active

  # --- Omega basis matrices ---------------------------------------------------
  Omega_sel_con <- matrix(rep(1, n), nrow = 1)
  Omega_sel_mod <- matrix(z, nrow = 1)

  Omega_out_con <- matrix(rep(1, n_active), nrow = 1)
  Omega_out_mod <- matrix(z[active], nrow = 1)

  Omega_out_con_full <- matrix(rep(1, n), nrow = 1)
  Omega_out_mod_full <- matrix(z, nrow = 1)

  # --- Prior SDs --------------------------------------------------------------
  Sigma0_sel_con <- matrix(sd_control_hurdle^2 / ntree_control, nrow = 1)
  Sigma0_sel_mod <- matrix(sd_moderate_hurdle^2 / ntree_moderate, nrow = 1)
  Sigma0_out_con <- matrix(con_sd_out^2 / ntree_control, nrow = 1)
  Sigma0_out_mod <- matrix(mod_sd_out^2 / ntree_moderate, nrow = 1)

  # --- Chipman-style outcome-variance prior, calibrated on the active sample -
  if (is.null(lambda)) {
    y_scaled_active <- y_cpp[active]
    lmf <- tryCatch(
      suppressWarnings(lm(y_scaled_active ~ x_c_out_active)),
      error = function(e) NULL
    )
    sighat_scaled <- if (!is.null(lmf)) {
      s <- summary(lmf)$sigma
      if (is.na(s) || s <= 0) 0.9 else s
    } else {
      0.9
    }
    qchi <- qchisq(1.0 - sigq, nu)
    lambda <- (sighat_scaled^2 * qchi) / nu
  }

  # --- Call C++ core ----------------------------------------------------------
  res <- pathca_bcfCore(y_ = y_cpp,
                        I_ = as.numeric(active),
                        Omega_sel_con = Omega_sel_con, Omega_sel_mod = Omega_sel_mod,
                        Omega_out_con = Omega_out_con, Omega_out_mod = Omega_out_mod,
                        Omega_out_con_full = Omega_out_con_full, Omega_out_mod_full = Omega_out_mod_full,
                        Omega_sel_con_est = Omega_sel_con_est, Omega_sel_mod_est = Omega_sel_mod_est,
                        Omega_out_con_est = Omega_out_con_est, Omega_out_mod_est = Omega_out_mod_est,
                        x_sel_con_ = t(x_c_sel), x_sel_mod_ = t(x_m_sel),
                        x_out_con_ = t(x_c_out_active), x_out_mod_ = t(x_m_out_active),
                        x_out_con_full_ = t(x_c_out_full), x_out_mod_full_ = t(x_m_out_full),
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

  # --- Re-scale outcome posterior draws back to log scale ---------------------
  res$out_con_post <- muy_active + sdy_active * res$out_con_post
  res$out_mod_post <- sdy_active * res$out_mod_post
  res$yhat_post <- muy_active + sdy_active * res$yhat_post

  res$beta_post <- sdy_active * res$beta_post
  res$sigma0_post <- sdy_active * res$sigma0_post
  res$sigma_post <- sdy_active * res$sigma_post

  res$sel_tau_post <- res$sel_tau_post
  res$out_tau_post <- sdy_active * res$out_tau_post

  res$V_post <- muy_active + sdy_active * res$V_post

  res$pihat_sel <- pihat_sel
  res$pihat_active <- pihat_active
  res$muy_active <- muy_active
  res$sdy_active <- sdy_active

  if (!is.null(x_control_est)) {
    res$out_con_est_post <- muy_active + sdy_active * res$out_con_est_post
    res$out_mod_est_post <- sdy_active * res$out_mod_est_post
    res$yhat_est_post <- res$out_con_est_post + res$out_mod_est_post
    res$pihat_sel_est <- pihat_sel_est
    res$pihat_active_est <- pihat_active_est
  }

  return(res)
}
