#   matrixt.R
#   MixMatrix: Classification with Matrix Variate Normal and t distributions
#   Copyright (C) 2018-9  GZ Thompson <gzthompson@gmail.com>
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#   along with this program; if not, a copy is available at
#   https://www.R-project.org/Licenses/


#' Distribution functions for the matrix variate t distribution.
#'
#' Density and random generation for the matrix variate t distribution.
#' @name rmatrixt
#' @param n number of observations for generation
#' @param x quantile for density
#' @param df  degrees of freedom (\eqn{>0}, may be non-integer),
#'    `df = 0, Inf` is allowed and will return a normal distribution.
#' @param mean \eqn{p \times q}{p * q} This is really a 'shift' rather than a
#'    mean, though the expected value will be equal to this if
#'    \eqn{df > 2}
#' @param L \eqn{p \times p}{p * p}  matrix specifying relations among the rows.
#'     By default, an identity matrix.
#' @param R \eqn{q \times q}{q * q}  matrix specifying relations among the
#'     columns. By default, an identity matrix.
#' @param U \eqn{LL^T}  - \eqn{p \times p}{p * p}   positive definite matrix for
#'    rows, computed from \eqn{L} if not specified.
#' @param V \eqn{R^T R}  - \eqn{q \times q}{q * q}  positive definite matrix for
#'    columns, computed from \eqn{R}  if not specified.
#' @param list Defaults to `FALSE` . If this is `TRUE` , then the
#'    output will be a list of matrices.
#' @param array If \eqn{n = 1}  and this is not specified and `list`  is
#'    `FALSE` , the function will return a matrix containing the one
#'    observation. If \eqn{n > 1} , should be the opposite of `list` .
#'    If `list`  is `TRUE` , this will be ignored.
#' @param force In `rmatrix`: if `TRUE`, will take the input of
#'    `R` directly - otherwise uses `V` and uses Cholesky
#'    decompositions. Useful for generating degenerate t-distributions.
#'    Will also override concerns about potentially singular matrices
#'    unless they are not, in fact, invertible.
#' @param log logical; in `dmatrixt`, if `TRUE`, probabilities
#'    `p` are given as `log(p)`.
#' @return `rmatrixt` returns either a list of \eqn{n}
#'    \eqn{p \times q}{p * q}  matrices or a
#' \eqn{p \times q \times n}{p * q * n}
#'    array.
#'
#'    `dmatrixt` returns the density at `x`.
#' @references  Gupta, Arjun K, and Daya K Nagar. 1999.
#' Matrix Variate Distributions.
#'     Vol. 104. CRC Press. ISBN:978-1584880462
#'
#' Dickey, James M. 1967. “Matricvariate Generalizations of the Multivariate t
#'        Distribution and the Inverted Multivariate t
#'        Distribution.” Ann. Math. Statist. 38 (2): 511–18.
#' \doi{10.1214/aoms/1177698967}
#' @details
#' The matrix \eqn{t}-distribution is parameterized slightly
#'  differently from the univariate and multivariate \eqn{t}-distributions
#'  - the variance is scaled by a factor of `1/df`.
#'  In this parameterization, the variance for a \eqn{1 \times 1}{1 * 1} matrix
#'  variate \eqn{t}-distributed random variable with identity variance matrices
#'  is \eqn{1/(df-2)} instead of \eqn{df/(df-2)}. A Central Limit Theorem
#'  for the matrix variate \eqn{T} is then that as `df` goes to
#'  infinity, \eqn{MVT(0, df, I_p, df*I_q)} converges to
#'  \eqn{MVN(0,I_p,I_q)}.
#'
#' @seealso [rmatrixnorm()],
#'     [rmatrixinvt()],[rt()] and
#' [stats::Distributions()].
#'
#' @export
#'
#' @examples
#' set.seed(20180202)
#' # random matrix with df = 10 and the given mean and L matrix
#' rmatrixt(
#'   n = 1, df = 10, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
#'   L = matrix(c(2, 1, 0, .1), nrow = 2), list = FALSE
#' )
#' # comparing 1-D distribution of t to matrix
#' summary(rt(n = 100, df = 10))
#' summary(rmatrixt(n = 100, df = 10, matrix(0)))
#' # demonstrating equivalence of 1x1 matrix t to usual t
#' set.seed(20180204)
#' x <- rmatrixt(n = 1, mean = matrix(0), df = 1)
#' dt(x, 1)
#' dmatrixt(x, df = 1)
rmatrixt <- function(n, df, mean,
                     L = diag(dim(as.matrix(mean))[1]),
                     R = diag(dim(as.matrix(mean))[2]), U = L %*% t(L),
                     V = t(R) %*% R,
                     list = FALSE,
                     array = NULL,
                     force = FALSE) {
  if (!(allnumeric_stop(n, df, mean, L, R, U, V))) {
    stop("Non-numeric input. ")
  }
  if (!(n > 0)) {
    stop("n must be > 0. n =", n)
  }

  if (length(df) != 1) stop("Length of df must be 1. length = ", length(df))

  if (((is.null(df)) || is.na(df) || (df < 0))) {
    stop("df must be >= 0. df =", df)
  }
  if (df == 0 || is.infinite(df)) {
    return(rmatrixnorm(
      n = n, mean = mean, U = U,
      V = V, list = list, array = array
    ))
  }
  mean <- as.matrix(mean)
  u_mat <- as.matrix(U)
  v_mat <- as.matrix(V)
  dims <- dim(mean)
  dim_result <- dimcheck_stop(u_mat, v_mat, dims)
  if (!(dim_result == "ok")) stop(dim_result)
  if (force && !missing(R)) chol_v <- R else chol_v <- chol(v_mat)

  if (min(diag(chol_v)) < 1e-6 && !force) {
    stop(
      "Potentially singular covariance, use force = TRUE if intended. ",
      min(diag(chol_v))
    )
  }

  nobs <- prod(dims) * n
  mat <- array(stats::rnorm(nobs), dim = c(dims, n))

  chol_u <- CholWishart::rInvCholWishart(n, df + dims[1] - 1, u_mat)

  result <- array(dim = c(dims, n))

  for (i in seq(n)) {
    result[, , i] <- mean + (crossprod(chol_u[, , i], mat[, , i])) %*% (chol_v)
  }

  if (n == 1 && list == FALSE && is.null(array)) {
    return(array(result[, , 1], dim = dims[1:2]))
    # if n = 1 and you don't specify arguments, it just returns a matrix
  }
  if (list) {
    return(lapply(seq(dim(result)[3]), function(x) result[, , x]))
  }
  if (!(list) && !(is.null(array))) {
    if (!(array)) {
      warning("list FALSE and array FALSE, returning array")
    }
  }
  return(result)
}

#' @rdname rmatrixt
#' @export
dmatrixt <- function(x, df, mean = matrix(0, p, n),
                     L = diag(p),
                     R = diag(n), U = L %*% t(L),
                     V = t(R) %*% R, log = FALSE) {
  dims <- dim(x)
  if (is.null(dims) || length(dims) == 1) x <- matrix(x)

  dims <- dim(x)
  if (length(dims) == 2) x <- array(x, dim = (dims <- c(dims, 1)))
  p <- dims[1]
  n <- dims[2]
  if (!(allnumeric_stop((x), df, (mean), (L), (R), (U), (V)))) {
    stop("Non-numeric input. ")
  }
  if (length(df) != 1) stop("Length of df must be 1. length = ", length(df))
  if (((is.null(df)) || is.na(df) || (df < 0))) {
    stop("df must be >= 0. df =", df)
  }

  mean <- as.matrix(mean)
  u_mat <- as.matrix(U)
  v_mat <- as.matrix(V)
  dim_result <- dimcheck_stop(u_mat, v_mat, dims)
  if (!(dim_result == "ok")) stop(dim_result)
  if ((df == 0 || is.infinite(df))) {
    return(dmatrixnorm(x, mean = mean, U = u_mat, V = v_mat, log = log))
  }

  # gammas is constant
  # this could be shifted into C++ but I don't want to pull out of CholWishart
  gammas <- as.numeric(CholWishart::lmvgamma(
    (0.5) * (df + dims[1] + dims[2] - 1),
    dims[1]
  ) -
    0.5 * dims[1] * dims[2] * log(pi) -
    CholWishart::lmvgamma(0.5 * (df + dims[1] - 1), dims[1]))

  results <- as.numeric(dmat_t_calc(x, df, mean, u_mat, v_mat))
  results <- results + gammas
  if (log) {
    return(results)
  } else {
    return(exp(results))
  }
}


#' @title Distribution functions for matrix variate inverted t distributions
#'
#' @description Generate random samples from the inverted matrix
#'    variate t distribution or compute densities.
#'
#' @name rmatrixinvt
#' @rdname rmatrixinvt
#' @inheritParams rmatrixt
#' @return `rmatrixinvt` returns either a list of \eqn{n}
#'    \eqn{p \times q}{p * q}  matrices or
#'    a \eqn{p \times q \times n}{p * q * n}  array.
#'
#'    `dmatrixinvt` returns the density at  `x`.
#'
#' @seealso  [rmatrixnorm()], [rmatrixt()],
#'    and [stats::Distributions()].
#'
#' @references Gupta, Arjun K, and Daya K Nagar. 1999.
#' Matrix Variate Distributions.
#'     Vol. 104. CRC Press. ISBN:978-1584880462
#'
#'   Dickey, James M. 1967. “Matricvariate Generalizations of the Multivariate t
#'        Distribution and the Inverted Multivariate t
#'        Distribution.” Ann. Math. Statist. 38 (2): 511–18.
#' \doi{10.1214/aoms/1177698967}

#' @export
#' @examples
#' # an example of drawing from the distribution and computing the density.
#' A <- rmatrixinvt(n = 2, df = 10, diag(4))
#' dmatrixinvt(A[, , 1], df = 10, mean = diag(4))
rmatrixinvt <- function(n, df, mean,
                        L = diag(dim(as.matrix(mean))[1]),
                        R = diag(dim(as.matrix(mean))[2]),
                        U = L %*% t(L), V = t(R) %*% R,
                        list = FALSE, array = NULL) {
  if (!(allnumeric_stop(n, df, (mean), (L), (R), (U), (V)))) {
    stop("Non-numeric input. ")
  }
  if (length(df) != 1) stop("Length of df must be 1. length = ", length(df))
  if (((is.null(df)) || is.na(df) || (df < 0))) {
    stop("df must be >= 0. df =", df)
  }
  if (!(n > 0)) {
    stop("n must be > 0.", n)
  }
  mean <- as.matrix(mean)
  u_mat <- as.matrix(U)
  v_mat <- as.matrix(V)

  dims <- dim(mean)
  # checks for conformable matrix dimensions
  dim_result <- dimcheck_stop(u_mat, v_mat, dims)
  if (!(dim_result == "ok")) stop(dim_result)
  nobs <- prod(dims) * n
  mat <- array(stats::rnorm(nobs), dim = c(dims, n))

  s_mat <- stats::rWishart(n, df + dims[1] - 1, diag(dims[1]))

  result <- rmat_inv_t_calc(s_mat, mat, u_mat, v_mat, mean)


  if (n == 1 && list == FALSE && is.null(array)) {
    return(array(result[, , 1], dim = dims[1:2]))
    # if n = 1 and you don't specify arguments, it just returns a matrix
  }
  if (list) {
    return(lapply(seq(dim(result)[3]), function(x) result[, , x]))
  }
  if (!(list) && !(is.null(array))) {
    if (!(array)) {
      warning("list FALSE and array FALSE, returning array")
    }
  }
  (result)
}



#' @rdname rmatrixinvt
#' @export
dmatrixinvt <- function(x, df, mean = matrix(0, p, n), L = diag(p),
                        R = diag(n), U = L %*% t(L), V = t(R) %*% R,
                        log = FALSE) {
  dims <- dim(x)
  if (is.null(dims) || length(dims) == 1) x <- matrix(x)

  dims <- dim(x)
  if (length(dims) == 2) x <- array(x, dim = (dims <- c(dims, 1)))
  p <- dims[1]
  n <- dims[2]
  if (!(allnumeric_stop((x), df, (mean), (L), (R), (U), (V)))) {
    stop("Non-numeric input. ")
  }
  if (length(df) != 1) stop("Length of df must be 1. length = ", length(df))
  if (((is.null(df)) || is.na(df) || (df < 0))) {
    stop("df must be >= 0. df =", df)
  }

  mean <- as.matrix(mean)
  u_mat <- as.matrix(U)
  v_mat <- as.matrix(V)
  dim_result <- dimcheck_stop(u_mat, v_mat, dims)
  if (!(dim_result == "ok")) stop(dim_result)
  gammas <- as.numeric(
    CholWishart::lmvgamma((0.5) * (df + dims[1] + dims[2] - 1), dims[1]) -
      0.5 * prod(dims[1:2]) * log(pi) -
      CholWishart::lmvgamma(0.5 * (df + dims[1] - 1), dims[1])
  )

  results <- as.numeric(dmat_inv_t_calc(x, df, mean, u_mat, v_mat))
  if (any(is.nan(results))) {
    warning("warning: probability distribution undefined when det < 0.")
  }
  results <- gammas + results
  if (log) {
    return(results)
  } else {
    return(exp(results))
  }
}
