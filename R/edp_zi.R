## EDP+ZI benchmark (Kim et al. 2024), adapted for the zicbcf simulation DGPs:
## binary treatment + FIVE continuous covariates. Pure-R driver around the
## compiled edp_clustering_c() sampler.

## --- helpers -------------------------------------------------------------
.edp_rinvchisq <- function(n, df, scale = 1/df) {
  df <- rep(df, len = n); scale <- rep(scale, len = n)
  (df * scale) / rchisq(n, df = df)
}
## expit/logistic link of the model, exp(eta)/(1+exp(eta)); evaluated with
## plogis() (== the same function) to avoid exp() overflow for large |eta|.
.edp_expit <- function(beta, x) plogis(as.vector(x %*% beta))
.edp_rdir <- function(alpha) { g <- rgamma(length(alpha), alpha, 1); g / sum(g) }
## E[X | X > 0] for X ~ N(m, s): the truncated-normal mean m + s*phi(a)/(1-Phi(a))
## with a = -m/s, evaluated via the log-space Mills ratio so the deep-tail
## 0/0 (which makes truncnorm::etruncnorm return NaN) is avoided. Same quantity.
.edp_etrunc <- function(m, s) {
  m <- as.vector(m)
  if (s <= 0) return(pmax(m, 0))
  a <- -m / s
  res <- m + s * exp(dnorm(a, log = TRUE) - pnorm(a, lower.tail = FALSE, log.p = TRUE))
  ## exact limits for near-degenerate clusters (tiny s -> |a| huge): the
  ## Mills ratio formula overflows/cancels, but the truncated mean is exactly
  ## 0 when m is far below the truncation point and exactly m when far above.
  hi <- a > 37; lo <- a < -37
  if (any(hi)) res[hi] <- 0
  if (any(lo)) res[lo] <- m[lo]
  res
}

#' Fit the EDP+ZI benchmark (Kim et al. 2024) for zero-inflated causal effects
#'
#' Enriched Dirichlet process mixture for heterogeneous treatment effects with
#' zero-inflated outcomes. This is a faithful re-implementation of the authors'
#' EDPcausal sampler (\code{clustering.cpp} + \code{function.R}) adapted to the
#' zicbcf simulation DGPs: the binary treatment uses a categorical kernel and
#' all five covariates use Gaussian kernels (the paper's prescribed handling of
#' continuous covariates). The MODEL (EDP mixture, kernels, priors, the
#' zero-inflation link \eqn{\pi=\mathrm{expit}(r)}, and the G-computation of
#' \eqn{E[Y\mid do(a),x]}) is exactly as published; implementation defects in the
#' reference code (a \code{1/e} constant in place of the expit, a frozen
#' concentration-parameter update, and product-space under/overflow) are
#' corrected without altering the model.
#'
#' @param y Semicontinuous response (>= 0).
#' @param z Binary treatment assignment (0/1).
#' @param x Numeric covariate matrix (n x 5).
#' @param nburn Burn-in iterations (default 1000).
#' @param nsim Post-burn-in iterations retained (default 1000).
#' @param nthin Thinning interval for the retained draws (default 1).
#' @param alpha_theta,alpha_omega Initial DP / nested-DP concentrations
#'   (Gamma(1,1) priors; defaults 2).
#' @param Ky_init Number of initial outcome clusters (default 5).
#' @param update_interval Iterations between progress messages (0 = silent).
#' @return A list with \code{cate} (draws x n matrix of
#'   \eqn{E[Y(1)-Y(0)\mid x]}), \code{ate} (row means), and \code{p1}, \code{p0}
#'   (draws x n matrices of \eqn{P(Y>0\mid do(a),x)} for a = 1, 0).
#' @export
edp_zi <- function(y, z, x,
                   nburn = 1000, nsim = 1000, nthin = 1,
                   alpha_theta = 2, alpha_omega = 2,
                   Ky_init = 5, update_interval = 200) {

  n <- length(y)
  Xmat <- as.matrix(x); storage.mode(Xmat) <- "double"
  a <- as.numeric(z)
  x1 <- Xmat[,1]; x2 <- Xmat[,2]; x3 <- Xmat[,3]; x4 <- Xmat[,4]; x5 <- Xmat[,5]
  xa <- cbind(1, a, x1, x2, x3, x4, x5); storage.mode(xa) <- "double"

  ## hyper-parameters (match paper / function.R)
  alpha00 <- c(0.1, 0.1)
  nu0 <- 2; tau0 <- 1; c0 <- 1; mu0 <- 0
  mu.r <- 0; Sigma.r <- 5; mu.beta <- 0; Sigma.beta <- 10

  ## ---- adapted update_parameters (5 continuous covariates) ----
  ## The concentration parameters are updated with the standard Escobar & West
  ## (1995) augmentation under the paper's Gamma(1,1) priors, using the CURRENT
  ## alpha values (threaded in). The authors' function.R instead reads fixed
  ## globals, which never advances the alpha chain -- corrected here without
  ## changing the prior or the augmentation scheme.
  update_parameters <- function(Sx, Sy, xPAR_p0, xPAR_mu, xPAR_sig,
                                yPAR_beta, yPAR_r, yPAR_sig,
                                alpha_theta, alpha_omega) {
    Sy <- as.vector(Sy); Sx <- as.vector(Sx)
    numY <- dim(yPAR_beta)[2]
    sysx <- cbind(Sy, Sx)
    uniqueYX <- unique(sysx)
    uniqueYX <- uniqueYX[order(uniqueYX[,1], uniqueYX[,2]), , drop = FALSE]
    uniqueY <- unique(Sy)

    for (j in 1:length(uniqueY)) {
      idxj  <- (Sy == j)
      X_lj  <- xa[idxj, , drop = FALSE]
      ysub  <- y[idxj]

      ## --- sigma^2 of the truncated-normal outcome (MH) ---
      e_r <- .edp_expit(yPAR_r[, j], X_lj)
      prop <- rgamma(1, 10000 * yPAR_sig[j, 1], 10000)
      num <- sum(log(e_r * (ysub == 0) + (ysub != 0) * (1 - e_r) *
                 dtruncnorm(ysub, a = 0, b = Inf, mean = X_lj %*% yPAR_beta[, j], sd = sqrt(prop)))) +
             dgamma(prop, 0.1, 10, log = TRUE) + dgamma(yPAR_sig[j, 1], 10000 * prop, 10000, log = TRUE)
      den <- sum(log(e_r * (ysub == 0) + (ysub != 0) * (1 - e_r) *
                 dtruncnorm(ysub, a = 0, b = Inf, mean = X_lj %*% yPAR_beta[, j], sd = sqrt(yPAR_sig[j, 1])))) +
             dgamma(yPAR_sig[j, 1], 0.1, 10, log = TRUE) + dgamma(prop, 10000 * yPAR_sig[j, 1], 10000, log = TRUE)
      if (log(runif(1)) < (num - den) & !is.nan(num - den)) yPAR_sig[j, 1] <- prop

      ## --- regression coefficients beta (MH, one at a time) ---
      for (cc in 1:length(yPAR_beta[, j])) {
        prop_b <- yPAR_beta[, j]; prop_b[cc] <- yPAR_beta[cc, j] + rnorm(1, 0, 0.05)
        num <- sum(log(e_r * (ysub == 0) + (ysub != 0) * (1 - e_r) *
                   dtruncnorm(ysub, a = 0, b = Inf, mean = X_lj %*% prop_b, sd = sqrt(yPAR_sig[j, 1])))) +
               dnorm(prop_b[cc], mu.beta, sqrt(Sigma.beta), log = TRUE)
        den <- sum(log(e_r * (ysub == 0) + (ysub != 0) * (1 - e_r) *
                   dtruncnorm(ysub, a = 0, b = Inf, mean = X_lj %*% yPAR_beta[, j], sd = sqrt(yPAR_sig[j, 1])))) +
               dnorm(yPAR_beta[cc, j], mu.beta, sqrt(Sigma.beta), log = TRUE)
        if (log(runif(1)) < (num - den) & !is.nan(num - den)) yPAR_beta[cc, j] <- prop_b[cc]
      }

      ## --- zero-inflation coefficients r (MH, one at a time) ---
      for (cc in 1:length(yPAR_r[, j])) {
        prop_r <- yPAR_r[, j]; prop_r[cc] <- yPAR_r[cc, j] + rnorm(1, 0, 0.05)
        e_rp <- .edp_expit(prop_r, X_lj)
        num <- sum(log(e_rp * (ysub == 0) + (ysub != 0) * (1 - e_rp) *
                   dtruncnorm(ysub, a = 0, b = Inf, mean = X_lj %*% yPAR_beta[, j], sd = sqrt(yPAR_sig[j, 1])))) +
               dnorm(prop_r[cc], mu.r, sqrt(Sigma.r), log = TRUE)
        den <- sum(log(e_r * (ysub == 0) + (ysub != 0) * (1 - e_r) *
                   dtruncnorm(ysub, a = 0, b = Inf, mean = X_lj %*% yPAR_beta[, j], sd = sqrt(yPAR_sig[j, 1])))) +
               dnorm(yPAR_r[cc, j], mu.r, sqrt(Sigma.r), log = TRUE)
        if (log(runif(1)) < (num - den) & !is.nan(num - den)) { yPAR_r[cc, j] <- prop_r[cc]; e_r <- e_rp }
      }

      ## --- covariate kernels within each X-subcluster of Y-cluster j ---
      nXj <- length(which(uniqueYX[,1] == j))
      for (l in 1:nXj) {
        ind  <- which(uniqueYX[,1] == j & uniqueYX[,2] == l)
        cell <- (Sy == j & Sx == l)
        ncell <- sum(cell)
        ## treatment kernel (binary, Dirichlet-multinomial)
        xPAR_p0[ind, ] <- .edp_rdir(alpha00 + as.numeric(table(factor(a[cell] + 1, levels = 1:2))))
        ## five continuous covariate kernels (Normal / scaled-inv-chisq)
        for (cc in 1:5) {
          xv  <- Xmat[cell, cc]
          sd2 <- ifelse(is.na(stats::sd(xv)^2), 0, stats::sd(xv)^2)
          xPAR_sig[cc, ind] <- .edp_rinvchisq(1, nu0 + ncell,
              (nu0 * tau0 + (ncell - 1) * sd2 + c0 * ncell / (c0 + ncell) * (mean(xv) - mu0)^2) / (nu0 + ncell))
          term1 <- c0 / tau0
          term2 <- ncell / xPAR_sig[cc, ind]
          xPAR_mu[cc, ind] <- rnorm(1, (term1 * mu0 + term2 * mean(xv)) / (term1 + term2), sqrt(1 / (term1 + term2)))
        }
      }
    }

    ## --- concentration parameters (Escobar & West augmentation) ---
    nj_i  <- table(factor(Sy, levels = 1:length(uniqueY)))
    nlj_i <- table(Sy, Sx)
    numXj <- apply(nlj_i, 1, function(v) length(which(v != 0)))

    eta <- rbeta(1, alpha_theta + 1, n)
    pival <- (1 + length(uniqueY) - 1) / (n * (1 - log(eta)))
    pival <- pival / (1 + pival)
    pi.temp <- sample(c(1, 0), 1, prob = c(pival, 1 - pival))
    alpha_theta <- pi.temp * rgamma(1, 1 + length(uniqueY), 1 - log(eta)) +
                   (1 - pi.temp) * rgamma(1, 1 + length(uniqueY) - 1, 1 - log(eta))

    prop <- rgamma(1, alpha_omega * 1000, 1000)
    rat <- log(dgamma(prop, 1, 1) * prop^(sum(numXj - 1)) * prod((prop + nj_i) * beta(prop + 1, nj_i))) +
           log(dgamma(alpha_omega, prop * 1000, 1000)) -
           log(dgamma(alpha_omega, 1, 1) * alpha_omega^(sum(numXj - 1)) * prod((alpha_omega + nj_i) * beta(alpha_omega + 1, nj_i))) -
           log(dgamma(prop, alpha_omega * 1000, 1000))
    if (log(runif(1)) <= rat & !is.nan(rat)) alpha_omega <- prop

    list(alpha_omega = alpha_omega, alpha_theta = alpha_theta,
         xPAR_p0 = xPAR_p0, xPAR_mu = xPAR_mu, xPAR_sig = xPAR_sig,
         yPAR_beta = yPAR_beta, yPAR_r = yPAR_r, yPAR_sig = yPAR_sig)
  }

  ## ---- initialisation ----
  Ky_init <- max(2, min(Ky_init, floor(n / 5)))
  Sy <- sample(rep(1:Ky_init, length.out = n))      # all labels guaranteed present
  Sx <- rep(1L, n)                                  # one X-subcluster per Y-cluster
  numY <- Ky_init; nX <- Ky_init
  state <- list(
    Sy = Sy, Sx = Sx,
    xPAR_p0  = matrix(0.5, nrow = nX, ncol = 2),
    xPAR_mu  = matrix(rnorm(5 * nX), nrow = 5, ncol = nX),
    xPAR_sig = matrix(1, nrow = 5, ncol = nX),
    yPAR_beta = matrix(rnorm(7 * numY), nrow = 7, ncol = numY),
    yPAR_r    = matrix(rnorm(7 * numY), nrow = 7, ncol = numY),
    yPAR_sig  = matrix(1, nrow = numY, ncol = 1),
    alpha_theta = alpha_theta, alpha_omega = alpha_omega)

  total <- nburn + nsim
  keep <- seq(nburn + 1, total, by = nthin)
  ndraw <- length(keep)
  cate_draws <- matrix(NA_real_, nrow = ndraw, ncol = n)
  p1_draws   <- matrix(NA_real_, nrow = ndraw, ncol = n)
  p0_draws   <- matrix(NA_real_, nrow = ndraw, ncol = n)
  di <- 0L

  X1 <- cbind(1, 1, x1, x2, x3, x4, x5)   # potential-outcome design under A=1
  X0 <- cbind(1, 0, x1, x2, x3, x4, x5)   # under A=0

  ## per-draw G-computation of E[Y|do(a),x] and P(Y=0|do(a),x)
  gcomp <- function(st, Xd, a_val) {
    Sy <- as.vector(st$Sy); Sx <- as.vector(st$Sx)
    numY <- ncol(st$yPAR_beta)
    nj_i  <- table(factor(Sy, levels = 1:numY))
    nlj_i <- table(factor(Sy, levels = 1:numY), Sx)
    numXj <- apply(nlj_i, 1, function(v) length(which(v != 0)))
    Nlj <- matrix(nrow = numY, ncol = 2)
    Nlj[1, ] <- c(1, numXj[1])
    if (numY != 1) for (cc in 2:numY) Nlj[cc, ] <- c(cumsum(numXj)[cc - 1] + 1, cumsum(numXj)[cc])
    nnlj <- as.numeric(t(nlj_i)); nnlj <- nnlj[nnlj != 0]
    pcol <- if (a_val == 1) 2L else 1L

    ## Membership probability prob[i, cy] is the authors' mixture weight
    ##   nj_i[cy]/(alpha_theta+n) * sum_x nnlj[x]/(alpha_omega+nj_i[cy]) *
    ##                              dcat(a) * prod_k dnorm(x_k | mu, sig),
    ## normalised across cy (their `prob <- prob/rowSums(prob)`). We evaluate
    ## the identical quantity in log-space (log-sum-exp) so that outlier
    ## subjects whose Gaussian kernels underflow to 0 do not yield 0/0 = NaN.
    log_un <- matrix(-Inf, nrow = n, ncol = numY)
    for (cy in 1:numY) {
      cols <- seq(Nlj[cy, 1], Nlj[cy, 2])
      Lx <- sapply(cols, function(xx)
        log(nnlj[xx] / (st$alpha_omega + nj_i[cy])) +
        log(st$xPAR_p0[xx, pcol]) +
        dnorm(Xd[,3], st$xPAR_mu[1, xx], sqrt(st$xPAR_sig[1, xx]), log = TRUE) +
        dnorm(Xd[,4], st$xPAR_mu[2, xx], sqrt(st$xPAR_sig[2, xx]), log = TRUE) +
        dnorm(Xd[,5], st$xPAR_mu[3, xx], sqrt(st$xPAR_sig[3, xx]), log = TRUE) +
        dnorm(Xd[,6], st$xPAR_mu[4, xx], sqrt(st$xPAR_sig[4, xx]), log = TRUE) +
        dnorm(Xd[,7], st$xPAR_mu[5, xx], sqrt(st$xPAR_sig[5, xx]), log = TRUE))
      if (is.null(dim(Lx))) Lx <- matrix(Lx, nrow = n)
      mx <- apply(Lx, 1, max)
      inner <- mx + log(rowSums(exp(Lx - mx)))                  # log-sum-exp over X-subclusters
      inner[!is.finite(mx)] <- -Inf                            # all X-subclusters impossible (e.g. no
                                                               # counterfactual-arm units) => weight 0
      log_un[, cy] <- log(nj_i[cy] / (st$alpha_theta + n)) + inner
    }
    mm <- apply(log_un, 1, max)
    bad <- !is.finite(mm)                                       # no admissible cluster (overlap edge) ->
    if (any(bad)) { log_un[bad, ] <- 0; mm[bad] <- 0 }          # fall back to uniform membership
    prob <- exp(log_un - mm)
    prob <- prob / rowSums(prob)                                # softmax across Y-clusters

    EYmat <- matrix(0, n, numY); PImat <- matrix(0, n, numY)
    for (cy in 1:numY) {
      Mu <- Xd %*% st$yPAR_beta[, cy]
      sig <- sqrt(st$yPAR_sig[cy, 1])
      pz <- .edp_expit(st$yPAR_r[, cy], Xd)        # P(structural zero)
      EYmat[, cy] <- .edp_etrunc(Mu, sig) * (1 - pz)
      PImat[, cy] <- pz
    }
    list(EY = rowSums(prob * EYmat), PI = rowSums(prob * PImat))
  }

  ## ---- MCMC ----
  for (l in 2:total) {
    clus <- edp_clustering_c(state$Sy, state$Sx, state$xPAR_p0, state$xPAR_mu, state$xPAR_sig,
                             state$yPAR_beta, state$yPAR_r, state$yPAR_sig,
                             state$alpha_omega, state$alpha_theta, xa, y)
    pars <- update_parameters(clus$Sx, clus$Sy, clus$xPAR_p0, clus$xPAR_mu, clus$xPAR_sig,
                              clus$yPAR_beta, clus$yPAR_r, clus$yPAR_sig,
                              state$alpha_theta, state$alpha_omega)
    state <- list(Sy = clus$Sy, Sx = clus$Sx,
                  xPAR_p0 = pars$xPAR_p0, xPAR_mu = pars$xPAR_mu, xPAR_sig = pars$xPAR_sig,
                  yPAR_beta = pars$yPAR_beta, yPAR_r = pars$yPAR_r, yPAR_sig = pars$yPAR_sig,
                  alpha_theta = pars$alpha_theta, alpha_omega = pars$alpha_omega)

    if (l > nburn && (l - nburn - 1) %% nthin == 0) {
      di <- di + 1L
      g1 <- gcomp(state, X1, 1)
      g0 <- gcomp(state, X0, 0)
      cate_draws[di, ] <- g1$EY - g0$EY
      p1_draws[di, ]   <- 1 - g1$PI
      p0_draws[di, ]   <- 1 - g0$PI
    }
    if (update_interval > 0 && l %% update_interval == 0)
      cat(sprintf("  [EDP] iter %d/%d  (#Y-clusters=%d)\n", l, total, ncol(state$yPAR_beta)))
  }

  cate_draws <- cate_draws[1:di, , drop = FALSE]
  p1_draws   <- p1_draws[1:di, , drop = FALSE]
  p0_draws   <- p0_draws[1:di, , drop = FALSE]
  list(cate = cate_draws, ate = rowMeans(cate_draws),
       p1 = p1_draws, p0 = p0_draws)
}
