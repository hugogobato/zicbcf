## Parametric Gamma benchmark models (Oganisian et al. 2019, appendix C):
## the "Gamma hurdle" two-part model and the naive "Gamma +.01" model, used
## there as parametric comparators for zero-inflated causal effects.
##
## The paper codes these in Stan with weakly-informative priors and reports
## posterior credible intervals. Here we use the standard Stan-free equivalent:
## a weakly-informative Bayesian GLM fit by penalized iteratively reweighted
## least squares (MAP estimate + Laplace / large-sample Gaussian posterior
## beta ~ N(beta_hat, V_hat)). The Gaussian prior plays the same role as the
## prior in the paper's Stan fit -- it keeps the naive +.01 model finite, whose
## unpenalized MLE diverges (a quasi-separation induced by the point mass that
## the +.01 trick maps to log(.01)). CATE / ATE / hurdle contrasts are obtained
## by G-computation (standardization), so the output objects are identical to
## the edp_zi() / dpglm() benchmarks:
##   cate -- draws x n matrix of E[Y(1)-Y(0) | x]
##   ate  -- row means of cate
##   p1, p0 -- draws x n matrices of P(Y > 0 | do(a), x), a = 1, 0
##             (NULL for Gamma +.01, which has no hurdle component).

## Build the fully treatment-interacted design (1, z, X, z*X) so the parametric
## model can express a covariate-varying treatment effect (needed for CATE).
.gamma_design <- function(zv, X) {
  Xi <- X * zv                          # column-wise z * covariate interactions
  colnames(Xi) <- paste0("zX", seq_len(ncol(X)))
  cbind(Intercept = 1, z = zv, X, Xi)
}

## symmetrise a covariance estimate (guards tiny numerical asymmetry that would
## otherwise trip MASS::mvrnorm).
.gamma_sym <- function(V) (V + t(V)) / 2

## Weakly-informative Bayesian GLM via penalized IRLS (Gaussian prior N(0, sd^2)
## on each coefficient). Returns the MAP coefficients and the Laplace posterior
## covariance phi * (X'WX + prior_prec)^{-1} (phi = 1 for binomial, Pearson
## dispersion for Gamma). The prior regularization makes the fit robust to the
## separation/divergence that defeats the plain IRLS MLE on the +.01 response.
.pen_irls <- function(D, y, family, prior_sd, maxit = 200, tol = 1e-8) {
  k <- ncol(D); n <- nrow(D)
  Pdiag <- 1 / prior_sd^2                 # diagonal prior precision (mean 0)
  beta <- rep(0, k)
  for (it in seq_len(maxit)) {
    eta <- pmin(pmax(as.vector(D %*% beta), -30), 30)   # clamp transient overflow
    mu  <- family$linkinv(eta)
    me  <- family$mu.eta(eta)
    v   <- pmax(family$variance(mu), 1e-10)
    w   <- as.vector(me^2 / v)
    zwk <- eta + (y - mu) / me            # IRLS working response
    A   <- crossprod(D, D * w); diag(A) <- diag(A) + Pdiag
    beta_new <- as.vector(solve(A, crossprod(D, w * zwk)))
    if (max(abs(beta_new - beta)) < tol * (1 + max(abs(beta)))) { beta <- beta_new; break }
    beta <- beta_new
  }
  eta <- pmin(pmax(as.vector(D %*% beta), -30), 30)
  mu  <- family$linkinv(eta)
  disp <- if (family$family == "Gamma")
    sum((y - mu)^2 / family$variance(mu)) / max(1, n - k) else 1   # Pearson dispersion
  A <- crossprod(D, D * as.vector(family$mu.eta(eta)^2 / pmax(family$variance(mu), 1e-10)))
  diag(A) <- diag(A) + Pdiag
  list(coef = beta, cov = .gamma_sym(disp * solve(A)))
}

## ndraws Gaussian-posterior draws of a fitted coefficient vector; always a matrix.
.gamma_post <- function(fit, ndraws)
  matrix(MASS::mvrnorm(ndraws, fit$coef, fit$cov), nrow = ndraws)

#' Fit the Gamma hurdle benchmark (Oganisian et al. 2019)
#'
#' Two-part parametric model: a logistic regression for the probability of a
#' positive outcome, \eqn{P(Y>0)}, and a Gamma regression with a log link for
#' the positive part \eqn{Y\mid Y>0}. Both parts use a fully treatment-by-covariate
#' interacted design \eqn{(1, z, X, z\!\cdot\!X)} so the model expresses a
#' covariate-varying treatment effect, and both are fit as weakly-informative
#' Bayesian GLMs (penalized IRLS: MAP estimate + Laplace posterior
#' \eqn{\beta\sim N(\hat\beta,\hat V)}), the Stan-free equivalent of the paper's
#' fit. Response-scale effects follow by G-computation,
#' \eqn{E[Y\mid do(a),x]=\mathrm{expit}(\eta_a)\cdot\exp(\xi_a)}.
#'
#' @param y Semicontinuous response (>= 0).
#' @param z Binary treatment assignment (0/1).
#' @param x Numeric covariate matrix (n x p).
#' @param nburn Ignored; kept for interface symmetry with the other benchmarks
#'   (the Laplace posterior needs no burn-in).
#' @param nsim Number of posterior draws to return (default 1000).
#' @param nthin Ignored; kept for interface symmetry.
#' @param prior_sd Prior standard deviation for the Gaussian coefficient prior
#'   (default 5, weakly informative).
#' @return A list with \code{cate} (draws x n matrix of \eqn{E[Y(1)-Y(0)\mid x]}),
#'   \code{ate} (row means), and \code{p1}, \code{p0} (draws x n matrices of
#'   \eqn{P(Y>0\mid do(a),x)} for a = 1, 0).
#' @export
gamma_hurdle <- function(y, z, x, nburn = 1000, nsim = 1000, nthin = 1, prior_sd = 5) {
  X <- as.matrix(x); storage.mode(X) <- "double"
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  n <- length(y); zv <- as.numeric(z); ndraws <- as.integer(nsim)

  D  <- .gamma_design(zv, X)
  D1 <- .gamma_design(rep(1, n), X)         # counterfactual designs for do(a)
  D0 <- .gamma_design(rep(0, n), X)

  ## part 1: logistic hurdle P(Y > 0) on the full sample
  fit_h <- .pen_irls(D, as.numeric(y > 0), stats::binomial(), prior_sd)
  ## part 2: Gamma(log) regression on the positive outcomes only
  pos   <- y > 0
  fit_g <- .pen_irls(D[pos, , drop = FALSE], y[pos], stats::Gamma(link = "log"), prior_sd)

  Bh <- .gamma_post(fit_h, ndraws)          # ndraws x k
  Bg <- .gamma_post(fit_g, ndraws)

  p1  <- stats::plogis(Bh %*% t(D1))        # ndraws x n
  p0  <- stats::plogis(Bh %*% t(D0))
  mu1 <- exp(pmin(Bg %*% t(D1), 30))
  mu0 <- exp(pmin(Bg %*% t(D0), 30))

  cate <- p1 * mu1 - p0 * mu0
  list(cate = cate, ate = rowMeans(cate), p1 = p1, p0 = p0)
}

#' Fit the naive Gamma +.01 benchmark (Oganisian et al. 2019)
#'
#' A single Gamma regression with a log link fit to the shifted outcome
#' \eqn{Y + 0.01} (the common, but ill-advised, trick of replacing structural
#' zeros with a small constant). The design is the fully treatment-by-covariate
#' interacted \eqn{(1, z, X, z\!\cdot\!X)}; the fit is a weakly-informative
#' Bayesian GLM (penalized IRLS) -- the prior is required here because the
#' unpenalized MLE diverges on the +.01 response. The response-scale CATE follows
#' by G-computation, \eqn{\tau(x)=\exp(\xi_1)-\exp(\xi_0)} (the +.01 cancels in
#' the difference). The model has no hurdle component, so it returns no
#' \eqn{P(Y>0)} contrast.
#'
#' @param y Semicontinuous response (>= 0).
#' @param z Binary treatment assignment (0/1).
#' @param x Numeric covariate matrix (n x p).
#' @param nburn Ignored; kept for interface symmetry with the other benchmarks.
#' @param nsim Number of posterior draws to return (default 1000).
#' @param nthin Ignored; kept for interface symmetry.
#' @param prior_sd Prior standard deviation for the Gaussian coefficient prior
#'   (default 5, weakly informative).
#' @param shift Constant added to the outcome (default 0.01).
#' @return A list with \code{cate} (draws x n matrix of \eqn{E[Y(1)-Y(0)\mid x]}),
#'   \code{ate} (row means), and \code{p1 = NULL}, \code{p0 = NULL} (no hurdle).
#' @export
gamma_plus01 <- function(y, z, x, nburn = 1000, nsim = 1000, nthin = 1,
                         prior_sd = 5, shift = 0.01) {
  X <- as.matrix(x); storage.mode(X) <- "double"
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  n <- length(y); zv <- as.numeric(z); ndraws <- as.integer(nsim)

  D  <- .gamma_design(zv, X)
  D1 <- .gamma_design(rep(1, n), X)
  D0 <- .gamma_design(rep(0, n), X)

  fit_g <- .pen_irls(D, y + shift, stats::Gamma(link = "log"), prior_sd)
  Bg <- .gamma_post(fit_g, ndraws)

  mu1 <- exp(pmin(Bg %*% t(D1), 30))        # E[Y + shift | do(1), x]
  mu0 <- exp(pmin(Bg %*% t(D0), 30))
  cate <- mu1 - mu0                         # shift cancels

  list(cate = cate, ate = rowMeans(cate), p1 = NULL, p0 = NULL)
}
