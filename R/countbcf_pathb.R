#' Fit a Unified Single-Forest Tweedie-BCF (Path B)
#'
#' \code{countbcf_pathb} fits the Tweedie Bayesian Causal Forest model (Path B)
#' for zero-inflated continuous / semicontinuous outcomes.
#' Under the hood, it uses an exact compound Poisson-Gamma (exponential) GIG-conjugate
#' formulation with Tweedie power index fixed at p = 1.5.
#'
#' @details
#' For a semicontinuous outcome \eqn{y_i \ge 0}, binary treatment indicator
#' \eqn{z_i \in \{0,1\}}, and covariates \eqn{x_i}:
#' \deqn{y_i \sim \mathrm{Tweedie}(\mu_i, \phi, p = 1.5)}
#' \deqn{\log \mu_i = \mu_f(x_i, \hat\pi_i) + z_i\, \tau_f(x_i),}
#' where \eqn{\mu_f} is a "vanilla" prognostic forest and \eqn{\tau_f} is a
#' moderating forest entered linearly in \eqn{z}.
#'
#' Leaf prior c and d parameters are calculated from a0 and a0_tau using the
#' standard GIG-prior hyperparameter mapping of Murray (2021).
#'
#' @param y Semicontinuous response variable (vector of non-negative doubles).
#' @param z Treatment assignment (vector of 0s and 1s).
#' @param x_control Design matrix for the prognostic function \eqn{\mu_f}.
#' @param x_moderate Design matrix for the treatment effect function \eqn{\tau_f}. Defaults to x_control.
#' @param pihat Length n estimate of \eqn{E(Z|X)}, included in mu-forest covariates. Defaults to 0.5.
#' @param offset Offset value (scalar or vector). If NULL, an offset of 1 is used.
#' @param nburn Number of burn-in MCMC iterations.
#' @param nsim Number of MCMC iterations to save after burn-in.
#' @param nthin Save every nthin'th MCMC iterate.
#' @param update_interval Print status every update_interval MCMC iterations.
#' @param ntree_control Number of trees in mu_f.
#' @param ntree_moderate Number of trees in tau_f.
#' @param a0 Concentration parameter for prognostic leaf prior. NA gives a default estimate.
#' @param a0_tau Concentration parameter for tau leaf prior. Default is a0/2.
#' @param base_control Base for tree prior on mu_f trees.
#' @param power_control Power for tree prior on mu_f trees.
#' @param base_moderate Base for tree prior on tau_f trees.
#' @param power_moderate Power for tree prior on tau_f trees.
#' @param include_pihat One of "control", "moderate", "both", "none".
#' @param randeff_design,randeff_variance_component_design,randeff_scales,randeff_df Random-effects setup; defaults disable random effects.
#' @param return_trees If TRUE, the serialized \code{tree_samples} for each
#'   forest are returned. Defaults to FALSE.
#' @param debug If TRUE, returns raw fit details for debugging.
#'
#' @return An object of class \code{"tweediebcf_fit"}: a list with posterior
#' draws of the log mean, per-forest fits, and dispersion parameter:
#' \describe{
#'   \item{\code{yhat_log}, \code{yhat}}{\eqn{n \times n_{\mathrm{sim}}} posterior
#'     draws of \eqn{\log E[Y_i \mid X_i, Z_i = z_i^{\mathrm{obs}}]} and its exponential.}
#'   \item{\code{order_vec}}{Permutation used to sort \code{y} (zeros first).
#'     Rows of \code{yhat*} follow this order.}
#'   \item{\code{mu_f_post}, \code{tau_f_post}}{
#'     \eqn{n_{\mathrm{sim}} \times n} per-iter raw coefficients
#'     of the prognostic / moderating forest in the \emph{original} input order.}
#'   \item{\code{phi}, \code{kappa}}{Posterior draws of the Tweedie dispersion parameter \eqn{\phi}.}
#' }
#'
#' @export
countbcf_pathb <- function(
  y,
  z,
  x_control,
  x_moderate = x_control,
  pihat = rep(0.5, length(y)),
  offset = NULL,
  nburn,
  nsim,
  nthin = 1,
  update_interval = 100,
  ntree_control = 250,
  ntree_moderate = 50,
  a0 = NA,
  a0_tau = NA,
  base_control = 0.95,
  power_control = 2,
  base_moderate = 0.25,
  power_moderate = 3,
  include_pihat = "control",
  randeff_design = matrix(1),
  randeff_variance_component_design = matrix(1),
  randeff_scales = 1,
  randeff_df = 3,
  return_trees = FALSE,
  debug = FALSE
) {

  # Fixed p=1.5 Tweedie model maps to model type 5 in countbcf_pathb C++ backend
  model_type <- 5

  if (any(is.na(y))) stop("Missing values in y")
  if (any(is.na(z))) stop("Missing values in z")
  if (any(is.na(x_control))) stop("Missing values in x_control")
  if (any(is.na(x_moderate))) stop("Missing values in x_moderate")
  if (any(is.na(pihat))) stop("Missing values in pihat")

  if (any(!is.finite(y))) stop("Non-numeric values in y")
  if (any(!is.finite(z))) stop("Non-numeric values in z")
  if (any(!is.finite(x_control))) stop("Non-numeric values in x_control")
  if (any(!is.finite(x_moderate))) stop("Non-numeric values in x_moderate")
  if (any(!is.finite(pihat))) stop("Non-numeric values in pihat")

  if (any(y < 0)) stop("Negative values in y")
  if (!all(sort(unique(z)) %in% c(0, 1))) stop("z must contain only 0 and 1")

  if (!identical(length(y), nrow(x_control)) || !identical(length(y), nrow(x_moderate)) ||
      !identical(length(y), length(z)) || !identical(length(y), length(pihat))) {
    stop("Data size mismatch among y, z, x_control, x_moderate, pihat")
  }

  if (nburn < 0) stop("nburn must be positive")
  if (nsim < 0) stop("nsim must be positive")
  if (nthin < 0) stop("nthin must be positive")
  if (nthin > nsim + 1) stop("nthin must be < nsim")
  if (nburn < 100) warning("A low (<100) value for nburn was supplied")

  pihat <- as.matrix(pihat)
  if (!(include_pihat %in% c("control", "moderate", "both", "none"))) {
    stop("include_pihat must be 'control', 'moderate', 'both', or 'none'")
  }
  if (include_pihat != "none" && length(unique(pihat)) == 1) {
    warning("All values of pihat are equal. pihat will not be included among covariates")
    include_pihat <- "none"
  }

  ### offset
  if (is.null(offset)) {
    mu0 <- rep(0, length(y))
  } else if (length(offset) == 1) {
    mu0 <- rep(log(offset), length(y))
  } else if (length(offset) == length(y)) {
    mu0 <- log(offset)
  } else {
    stop("offset incorrectly specified")
  }

  ### sort by y so the zeros come first (required by countbcf C++ code)
  order_vec <- order(y)
  y_sort <- y[order_vec]
  z_sort <- z[order_vec]
  mu0_sort <- mu0[order_vec]
  pihat_sort <- pihat[order_vec, , drop = FALSE]

  X_c <- matrix(x_control, ncol = ncol(x_control))[order_vec, , drop = FALSE]
  X_m <- matrix(x_moderate, ncol = ncol(x_moderate))[order_vec, , drop = FALSE]

  if (include_pihat %in% c("control", "both")) {
    X_c <- cbind(X_c, pihat_sort)
  }
  if (include_pihat %in% c("moderate", "both")) {
    X_m <- cbind(X_m, pihat_sort)
  }

  n <- length(y_sort)

  ### leaf hyperparameters
  if (is.na(a0)) {
    y_star <- quantile(y_sort, probs = 0.95)
    a0 <- 0.5 * (log(max(y_star, 1)) - mean(mu0_sort))
    if (!is.finite(a0) || a0 <= 0) a0 <- 1
  }
  if (is.na(a0_tau)) a0_tau <- a0 / 2

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

  z_c_mu <- 0; z_d_mu <- 0
  z_c_tau <- 0; z_d_tau <- 0

  ### leaf-scale (Sigma0) is unused in the loglinear leaf draw but required by the spec.
  Sigma0_dummy <- matrix(1, nrow = 1, ncol = 1)

  ### designs and specs
  X_list <- list(X_c, X_m)
  basis_list <- list(matrix(rep(1, n), ncol = 1), matrix(z_sort, ncol = 1))

  designs <- countbcf:::make_bart_designs(X_list, basis_list)

  make_spec <- function(design_idx, ntree, alpha, beta, vanilla, function_group) {
    design <- designs[[design_idx]]
    spec <- countbcf:::make_bart_spec(
      design = design,
      ntree = ntree,
      Sigma0 = Sigma0_dummy,
      scale_df = -1,
      vanilla = vanilla,
      alpha = alpha,
      beta = beta,
      update_leaf_scale = FALSE
    )
    spec$function_group <- function_group
    spec
  }

  specs <- list(
    make_spec(1, ntree_control,  base_control,  power_control,  TRUE,  0),  # mu_f
    make_spec(2, ntree_moderate, base_moderate, power_moderate, FALSE, 0)   # tau_f
  )

  for (k in seq_along(specs)) {
    is_tau <- !specs[[k]]$vanilla
    specs[[k]]$leaf_c <- if (is_tau) leaf_c_tau else leaf_c_mu
    specs[[k]]$leaf_d <- if (is_tau) leaf_d_tau else leaf_d_mu
  }

  lambda <- 1
  nu <- 3

  const_args <- list(
    y_     = y_sort,
    offset_ = mu0_sort,
    bart_specs = specs,
    bart_designs = designs,
    random_des = randeff_design,
    random_var_ix = randeff_variance_component_design,
    random_var = matrix(rep(1e-8, ncol(randeff_variance_component_design)), ncol = 1),
    random_var_df = randeff_df,
    randeff_scales = randeff_scales,
    burn = nburn, nd = nsim, thin = nthin,
    count_model = model_type,
    lambda = lambda, nu = nu,
    kappa_a = 1, kappa_b = 1, # weakly informative beta prior bounds
    leaf_c = leaf_c_mu, leaf_d = leaf_d_mu,
    z_c = z_c_mu, z_d = z_d_mu,
    kappa_prop_sd = 0.2,
    status_interval = update_interval,
    text_trace = TRUE
  )

  # Call the isolated C++ function
  fit_raw <- with(const_args, .Call(
    "_countbcf_countbcf_pathb",
    PACKAGE = "countbcf",
    y_, offset_, bart_specs, bart_designs,
    random_des, random_var, random_var_ix, random_var_df, randeff_scales,
    burn, nd, thin,
    count_model, lambda, nu,
    kappa_a, kappa_b, leaf_c, leaf_d, z_c, z_d,
    kappa_prop_sd,
    as.logical(return_trees),  # return_trees
    FALSE,                # save_trees
    FALSE,                # est_mod_fits
    FALSE,                # est_con_fits
    FALSE,                # prior_sample
    as.integer(status_interval),
    as.numeric(c(0.0)),   # lower_bd
    as.numeric(c(0.0)),   # upper_bd
    FALSE,                # probit
    text_trace,
    FALSE                 # R_trace
  ))

  ### organize per-forest fits
  forest_fits <- fit_raw$forest_fits

  out <- list(
    yhat_log = fit_raw$yhat_post,
    yhat = exp(fit_raw$yhat_post),
    order_vec = order_vec,
    sigma = fit_raw$sigma,
    phi = fit_raw$kappa,
    kappa = fit_raw$kappa
  )

  inv <- order(order_vec)
  n_sorted <- length(order_vec)
  .extract_coef <- function(cc) {
    nsim <- length(cc) / n_sorted
    arr <- matrix(as.numeric(cc), nrow = n_sorted, ncol = nsim)
    t(arr[inv, , drop = FALSE])         # nsim x n (original order)
  }

  out$mu_f_log  <- forest_fits[[1]]
  out$tau_f_log <- forest_fits[[2]]
  out$mu_f_post  <- .extract_coef(fit_raw$coefs[[1]])
  out$tau_f_post <- .extract_coef(fit_raw$coefs[[2]])

  if (return_trees) {
    out$control_fit  <- list(tree_samples = fit_raw$tree_trace[[1]])
    out$moderate_fit <- list(tree_samples = fit_raw$tree_trace[[2]])
  }

  if (debug) {
    out$const_args <- const_args
    out$raw_fit <- fit_raw
  }

  class(out) <- "tweediebcf_fit"
  return(out)
}
