#' zicbcf: Zero-Inflated Continuous Bayesian Causal Forests
#'
#' Zero-Inflated Continuous Bayesian Causal Forests (ZIC-BCF) and
#' ZIC-BCF with Duan's Smearing (ZIC-BCF-Smear) for semicontinuous outcomes.
#'
#' @docType package
#' @import Rcpp
#' @importFrom stats approxfun lm qchisq quantile sd
#' @importFrom Rcpp evalCpp
#' @useDynLib zicbcf, .registration = TRUE
#' @name zicbcf-package
#' @export tree_samples
NULL

loadModule("tree_samples", TRUE)
