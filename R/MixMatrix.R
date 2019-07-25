#' Functions for Matrix Distributions
#'
#' Provides sampling and density functions for matrix
#' variate normal, t, and inverted t distributions;  ML estimation for matrix
#' variate normal and t distributions including some restrictions on the
#' parameters using the EM algorithm; and classification by linear and quadratic discriminant
#' analysis for matrix variate normal and t distributions described
#' in \href{https://arxiv.org/abs/1907.09565}{Thompson et al. (2019)}.
#' Performs clustering with matrix variate normal and t mixture models.
#'
#' @docType package
#' @name MixMatrix
NULL
## usethis namespace: start
#' @useDynLib MixMatrix, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

