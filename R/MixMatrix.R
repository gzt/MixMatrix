#' Classification with Matrix Variate Normal and t Distributions
#'
#' Provides sampling and density functions for matrix
#' variate normal, \eqn{t}, and inverted \eqn{t} distributions;
#' ML estimation for matrix variate normal and \eqn{t} distributions
#' using the EM algorithm, including some restrictions on the parameters;
#' and classification by linear and quadratic discriminant
#' analysis for matrix variate normal and t distributions described
#' in [Thompson et al. (2019)](https://arxiv.org/abs/1907.09565).
#' Performs clustering with matrix variate normal and t mixture models.
#'
#' @keywords internal 
"_PACKAGE"

## usethis namespace: start
#' @useDynLib MixMatrix, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom glue glue_collapse
## usethis namespace: end
NULL
