#' Fit the Continuous Intensity Stage of Path D (Gamma BCF via GIG Conjugacy)
#'
#' Estimates the treatment-aware log-linear Gamma BCF model for positive outcomes
#' utilizing Subpopulation Propensity Adjustment (SPA) to control for selection-induced RIC bias.
#'
#' @param y Response variable (full vector or active subset of positive continuous outcomes).
#' @param z Treatment assignment (binary vector).
#' @param x_control Design matrix for prognostic covariates.
#' @param x_moderate Design matrix for moderating covariates.
#' @param pihat_pos Propensity score for the active subpopulation (SPA). If NULL, estimated via logistic regression on the active subset.
#' @param offset Offset vector.
#' @param nburn Number of burn-in MCMC iterations.
#' @param nsim Number of MCMC iterations to save.
#' @param nthin Thinning interval.
#' @param update_interval Print status every update_interval iterations.
#' @param ntree_control Number of trees in prognostic forest.
#' @param ntree_moderate Number of trees in moderating forest.
#' @param a0 Concentration parameter for prognostic leaf prior.
#' @param a0_tau Concentration parameter for moderating leaf prior.
#' @param base_control Base for tree prior on prognostic trees.
#' @param power_control Power for tree prior on prognostic trees.
#' @param base_moderate Base for tree prior on moderating trees.
#' @param power_moderate Power for tree prior on moderating trees.
#' @param kappa_a Prior shape alpha for shape parameter kappa_c.
#' @param kappa_b Prior shape beta for shape parameter kappa_c.
#' @param kappa_prop_sd SD of kappa M-H proposal.
#' @param include_pihat Whether to include pihat in "control", "moderate", "both", or "none".
#' @param return_trees Whether to return tree trace.
#'
#' @importFrom stats binomial glm optim predict sd
#' @export
pathd_gammabcf <- function(
  y,
  z,
  x_control,
  x_moderate = x_control,
  pihat_pos = NULL,
  offset = NULL,
  nburn = 500,
  nsim = 500,
  thin = 1,
  update_interval = 100,
  ntree_control = 250,
  ntree_moderate = 50,
  a0 = 1.0,
  a0_tau = 0.5,
  base_control = 0.95,
  power_control = 2,
  base_moderate = 0.25,
  power_moderate = 3,
  kappa_a = 5,
  kappa_b = 3,
  kappa_prop_sd = 0.15,
  include_pihat = "control",
  return_trees = FALSE
) {
  # 1. Filter to active subset where Y > 0
  active_idx <- which(y > 0)
  if (length(active_idx) == 0) {
    stop("No positive outcomes found in y.")
  }

  y_pos <- y[active_idx]
  z_pos <- z[active_idx]
  x_control_pos <- as.matrix(x_control)[active_idx, , drop = FALSE]
  x_moderate_pos <- as.matrix(x_moderate)[active_idx, , drop = FALSE]

  n_pos <- length(y_pos)

  # 2. SPA Propensity score adjustment
  if (is.null(pihat_pos)) {
    # Fit a simple logistic regression of Z on x_control within the active subset
    sd_cols <- apply(x_control_pos, 2, sd)
    keep_cols <- which(sd_cols > 1e-8)
    if (length(keep_cols) > 0) {
      fit_prop <- glm(z_pos ~ x_control_pos[, keep_cols, drop=FALSE], family = binomial(link = "logit"))
      pihat_pos <- predict(fit_prop, type = "response")
    } else {
      pihat_pos <- rep(mean(z_pos), n_pos)
    }
  } else {
    if (length(pihat_pos) == length(y)) {
      pihat_pos <- pihat_pos[active_idx]
    }
    if (length(pihat_pos) != n_pos) {
      stop("pihat_pos must have the same length as the active subset or the full dataset.")
    }
  }

  # 3. Offsets
  # The user supplies an exposure-style multiplicative offset c_i such that
  # log lambda_i = forest(x_i) + log c_i. The C++ backend stores
  # theta_i = -log lambda_i = offset_C + sum_trees, so offset_C = -log c_i.
  if (is.null(offset)) {
    mu0_pos <- rep(0, n_pos)
  } else {
    if (length(offset) == length(y)) {
      mu0_pos <- -log(offset[active_idx])
    } else if (length(offset) == 1) {
      mu0_pos <- rep(-log(offset), n_pos)
    } else if (length(offset) == n_pos) {
      mu0_pos <- -log(offset)
    } else {
      stop("Invalid offset length.")
    }
  }

  # 4. Include Propensity score in designs
  X_c <- x_control_pos
  X_m <- x_moderate_pos

  if (include_pihat %in% c("control", "both") && length(unique(pihat_pos)) > 1) {
    X_c <- cbind(X_c, pihat_pos)
  }
  if (include_pihat %in% c("moderate", "both") && length(unique(pihat_pos)) > 1) {
    X_m <- cbind(X_m, pihat_pos)
  }

  # 5. GIG Prior c, d parameters
  log_leaf_var_mu <- a0^2 / ntree_control
  c_start_mu <- ntree_control / a0^2 + 0.5
  leaf_c_mu <- optim(c_start_mu, fn = function(x) (trigamma(x) - log_leaf_var_mu)^2,
                     method = "Brent", lower = 0, upper = 1e9)$par
  leaf_d_mu <- exp(digamma(leaf_c_mu))

  log_leaf_var_tau <- a0_tau^2 / ntree_moderate
  c_start_tau <- ntree_moderate / a0_tau^2 + 0.5
  leaf_c_tau <- optim(c_start_tau, fn = function(x) (trigamma(x) - log_leaf_var_tau)^2,
                      method = "Brent", lower = 0, upper = 1e9)$par
  leaf_d_tau <- exp(digamma(leaf_c_tau))

  # 6. Designs and specs
  X_list <- list(X_c, X_m)
  basis_list <- list(matrix(rep(1, n_pos), ncol = 1), matrix(z_pos, ncol = 1))
  designs <- countbcf:::make_bart_designs(X_list, basis_list)

  Sigma0_dummy <- matrix(1, nrow = 1, ncol = 1)
  make_spec <- function(design_idx, ntree, alpha, beta, vanilla) {
    countbcf:::make_bart_spec(
      design = designs[[design_idx]],
      ntree = ntree,
      Sigma0 = Sigma0_dummy,
      scale_df = -1,
      vanilla = vanilla,
      alpha = alpha,
      beta = beta,
      update_leaf_scale = FALSE
    )
  }

  specs <- list(
    make_spec(1, ntree_control,  base_control,  power_control,  TRUE),  # mu_f
    make_spec(2, ntree_moderate, base_moderate, power_moderate, FALSE)  # tau_f
  )

  # Override priors via spec
  specs[[1]]$leaf_c <- leaf_c_mu
  specs[[1]]$leaf_d <- leaf_d_mu
  specs[[2]]$leaf_c <- leaf_c_tau
  specs[[2]]$leaf_d <- leaf_d_tau

  # 7. Call compiled backend
  fit_raw <- .Call(
    "_countbcf_pathd_gammabcf",
    PACKAGE = "countbcf",
    y_pos,
    mu0_pos,
    specs,
    designs,
    matrix(1),
    matrix(1e-8),
    matrix(1),
    3,
    1,
    nburn,
    nsim,
    thin,
    kappa_a,
    kappa_b,
    leaf_c_mu,
    leaf_d_mu,
    kappa_prop_sd,
    as.logical(return_trees),
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    as.integer(update_interval),
    TRUE
  )

  # 8. Organize and return results
  # The C++ backend stores theta = -log(lambda) for GIG conjugacy. Convert all
  # outputs to the conventional log-lambda scale so users can treat mu_f_post
  # and tau_f_post as the prognostic / treatment effects on log(lambda).
  log_lambda_post <- -fit_raw$yhat_post                       # nsim x n_pos
  lambda_post     <- exp(log_lambda_post)
  mu_f_post  <- -t(matrix(as.numeric(fit_raw$coefs[[1]]),
                          nrow = n_pos, ncol = nsim))
  tau_f_post <- -t(matrix(as.numeric(fit_raw$coefs[[2]]),
                          nrow = n_pos, ncol = nsim))

  lambda_0_post <- exp(mu_f_post)
  lambda_1_post <- exp(mu_f_post + tau_f_post)

  out <- list(
    yhat = lambda_post,
    yhat_log = log_lambda_post,
    mu_f_post = mu_f_post,
    tau_f_post = tau_f_post,
    lambda_0_post = lambda_0_post,
    lambda_1_post = lambda_1_post,
    tau_intensity_post = lambda_1_post - lambda_0_post,
    kappa = fit_raw$kappa,
    kappa_acpt = fit_raw$kappa_acceptance,
    active_idx = active_idx,
    pihat_pos = pihat_pos
  )

  if (return_trees) {
    out$control_fit <- list(tree_samples = fit_raw$tree_trace[[1]])
    out$moderate_fit <- list(tree_samples = fit_raw$tree_trace[[2]])
  }

  class(out) <- "pathd_gammabcf_fit"
  return(out)
}
