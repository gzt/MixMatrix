#   matrixnorm.R
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


#' Matrix variate Normal distribution functions
#'
#' @description  Density and random generation for the matrix variate
#' normal distribution
#'
#' @name rmatrixnorm
#' @param n number of observations to generate - must be a positive integer.
#' @param x quantile for density
#' @param mean \eqn{p \times q}{p * q}  matrix of means
#' @param L  \eqn{p \times p}{p * p}  matrix specifying relations
#' among the rows.
#'    By default, an identity matrix.
#' @param R \eqn{q \times q}{q * q}  matrix specifying relations among the
#' columns.  By default, an identity matrix.
#' @param U \eqn{LL^T}  - \eqn{p \times p}{p * p}
#' positive definite variance-covariance
#'    matrix for rows, computed from \eqn{L} if not specified.
#' @param V \eqn{R^T R}  - \eqn{q \times q}{q * q}
#' positive definite variance-covariance
#'    matrix for columns, computed from \eqn{R}  if not specified.
#' @param list Defaults to \code{FALSE} . If this is \code{TRUE} , then the
#'    output will be a list of matrices.
#' @param array If \eqn{n = 1}  and this is not specified and \code{list}  is
#'    \code{FALSE} , the function will return a matrix containing the one
#'    observation. If \eqn{n > 1} , should be the opposite of \code{list} .
#'    If \code{list}  is \code{TRUE} , this will be ignored.
#' @param force If TRUE, will take the input of \code{L} and/or \code{R}
#'    directly - otherwise computes \code{U} and \code{V} and uses Cholesky
#'    decompositions. Useful for generating degenerate normal distributions.
#'    Will also override concerns about potentially singular matrices
#'    unless they are not, in fact, invertible.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @return \code{rmatrixnorm} returns either a list of
#' \eqn{n}  \eqn{p \times q}{p * q}  matrices or
#'    a \eqn{p \times q \times n}{p * q * n}  array.
#'
#'    \code{dmatrixnorm} returns the density at \code{x}.
#'
#' @references Gupta, Arjun K, and Daya K Nagar. 1999.
#' Matrix Variate Distributions.
#'     Vol. 104. CRC Press. ISBN:978-1584880462

#' @seealso \code{\link{rmatrixt}}, \code{\link{rmatrixinvt}},
#'     \code{\link{rnorm}} and \code{\link[stats]{Distributions}}
#' @export
#' @examples
#' set.seed(20180202)
#' # a draw from a matrix variate normal with a certain mean
#' # and row-wise covariance
#' rmatrixnorm(
#'   n = 1, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
#'   L = matrix(c(2, 1, 0, .1), nrow = 2), list = FALSE
#' )
#' set.seed(20180202)
#' # another way of specifying this - note the output is equivalent
#' A <- rmatrixnorm(
#'   n = 10, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
#'   L = matrix(c(2, 1, 0, .1), nrow = 2), list = TRUE
#' )
#' A[[1]]
#' # demonstrating the dmatrixnorm function
#' dmatrixnorm(A[[1]],
#'   mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
#'   L = matrix(c(2, 1, 0, .1), nrow = 2), log = TRUE
#' )
rmatrixnorm <- function(n, mean,
                        L = diag(dim(as.matrix(mean))[1]),
                        R = diag(dim(as.matrix(mean))[2]),
                        U = L %*% t(L),
                        V = t(R) %*% R,
                        list = FALSE,
                        array = NULL,
                        force = FALSE) {
  if (!(all(
    is.numeric(n), is.numeric(mean), is.numeric(L), is.numeric(R),
    is.numeric(U), is.numeric(V)
  ))) {
    stop("Non-numeric input. ")
  }
  if (!(n > 0)) stop("n must be > 0. n = ", n)
  mean <- as.matrix(mean)
  U <- as.matrix(U)
  V <- as.matrix(V)
  if (missing(L)) {
    if (!symm.check(U)) stop("U not symmetric.")
  }
  if (missing(R)) {
    if (!symm.check(V)) stop("V not symmetric.")
  }
  dims <- dim(mean)

  # checks for conformable matrix dimensions
  if (!(dims[1] == dim(U)[2] && dim(U)[1] == dim(U)[2] &&
    dims[2] == dim(V)[1] && dim(V)[1] == dim(V)[2])) {
    stop("Non-conforming dimensions.", dims, dim(U), dim(V))
  }
  if (force && !missing(L)) cholU <- L else cholU <- chol.default(U)
  if (force && !missing(R)) cholV <- R else cholV <- chol.default(V)

  if (!force && (any(diag(cholU) < 1e-6) || any(diag(cholV) < 1e-6))) {
    stop(
      "Potentially singular covariance, use force = TRUE if intended. ",
      min(diag(cholU)), min(diag(cholV))
    )
  }
  ## this is the point at which C++ would take over
  nobs <- prod(dims) * n
  mat <- array(stats::rnorm(nobs), dim = c(dims, n))

  result <- array(apply(mat, 3, function(x) {
    mean +
      crossprod(cholU, x) %*% (cholV)
  }),
  dim = c(dims, n)
  )
  if (n == 1 && list == FALSE && is.null(array)) {
    return(result[, , 1])
    # if n = 1 and you don't specify arguments, if just returns a matrix
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


#' @rdname rmatrixnorm
#' @export
dmatrixnorm <- function(x, mean = matrix(0, p, n),
                        L = diag(p),
                        R = diag(n), U = L %*% t(L),
                        V = t(R) %*% R, log = FALSE) {
  # x <- as.matrix(x)
  dims <- dim(x)
  if (is.null(dims) || length(dims) == 1) x <- matrix(x)
  dims <- dim(x)
  if (length(dims) == 2) x <- array(x, dim = (dims <- c(dims, 1)))
  p <- dims[1]
  n <- dims[2]
  if (!(all(
    is.numeric(x), is.numeric(mean), is.numeric(L), is.numeric(R),
    is.numeric(U), is.numeric(V)
  ))) {
    stop("Non-numeric input. ")
  }

  mean <- as.matrix(mean)
  U <- as.matrix(U)
  V <- as.matrix(V)
  if (!symm.check(U)) stop("U not symmetric.")
  if (!symm.check(V)) stop("V not symmetric.")

  if (!(dims[1] == dim(U)[2] && dim(U)[1] == dim(U)[2] &&
    dims[2] == dim(V)[1] && dim(V)[1] == dim(V)[2])) {
    stop("Non-conforming dimensions.", dims, dim(U), dim(V))
  }

  logresult <- as.numeric(dmatnorm_calc(x, mean, U, V))

  if (log) {
    return(logresult)
  } else {
    return(exp(logresult))
  }
}

#' dmatrixnorm.unroll
#' @description Equivalent to dmatrixnorm except it works by unrolling
#'     to a vector. Alternatively, it can work on a matrix that has
#'     already been unrolled in the default R method (using
#'     \code{as.vector}), as data may be stored in that fashion.
#'
#'
#' @inheritParams rmatrixnorm
##'
#' @return Returns the density at the provided observation. This is an
#'    alternative method of computing which works by flattening out into
#'    a vector instead of a matrix.
#' @noRd
#' @keywords internal
#'
#'
#' @examples
#' set.seed(20180202)
#' A <- rmatrixnorm(
#'   n = 1, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
#'   L = matrix(c(2, 1, 0, .1), nrow = 2)
#' )
#' \dontrun{
#' dmatrixnorm.unroll(A,
#'   mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
#'   L = matrix(c(2, 1, 0, .1), nrow = 2), log = TRUE
#' )
#' }
#'
dmatrixnorm.unroll <- function(x, mean = array(0L, dim(as.matrix(x))),
                               L = diag(dim(mean)[1]), R = diag(dim(mean)[2]),
                               U = L %*% t(L), V = t(R) %*% R, log = FALSE,
                               unrolled = FALSE) {
  # results should equal other option - works by unrolling into MVN
  if (!(all(
    is.numeric(x), is.numeric(mean), is.numeric(L), is.numeric(R),
    is.numeric(U), is.numeric(V)
  ))) {
    stop("Non-numeric input. ")
  }
  x <- as.matrix(x)
  mean <- as.matrix(mean)
  U <- as.matrix(U)
  V <- as.matrix(V)
  if (!symm.check(U)) stop("U not symmetric.")
  if (!symm.check(V)) stop("V not symmetric.")
  dims <- dim(mean)
  if (unrolled) {
    dims <- c(dim(U)[1], dim(V)[1])
  }
  if (!(dims[1] == dim(U)[2] && dim(U)[1] == dim(U)[2] &&
    dims[2] == dim(V)[1] && dim(V)[1] == dim(V)[2])) {
    stop("Non-conforming dimensions. ", dims, dim(U), dim(V))
  }
  if (!unrolled) {
    vecx <- as.vector(x)
  } else {
    vecx <- x
  }
  meanx <- as.vector(mean)
  VU <- V %x% U
  # you should be using small enough matrices that determinants aren't a pain
  # also presumes not using a singular matrix normal distribution
  p <- dim(U)[1] # square matrices so only need first dimension
  n <- dim(V)[1]
  cholVU <- chol.default(VU)
  detVU <- prod(diag(cholVU))^2
  if (!(detVU > 1e-8)) stop("non-invertible matrix")
  UVinv <- chol2inv(cholVU)
  XM <- vecx - meanx
  logresult <- -0.5 * n * p * log(2 * pi) - 0.5 * log(detVU) -
    0.5 * sum(diag(crossprod(XM, UVinv) %*% XM))
  if (log) {
    return(logresult)
  } else {
    return(exp(logresult))
  }
}

#'  Maximum likelihood estimation for matrix normal distributions
#'
#' @description
#' Maximum likelihood estimates exist for \eqn{N > max(p/q,q/p)+1} and are
#' unique for \eqn{N > max(p,q)}. This finds the estimate for the mean and then
#' alternates between estimates for the \eqn{U} and \eqn{V} matrices until
#' convergence. An AR(1), compound symmetry, correlation matrix, or independence
#' restriction can be proposed for either or both variance matrices. However, if
#' they are inappropriate for the data, they may fail with a warning.
#'
#' @param data Either a list of matrices or a 3-D array with matrices in
#'    dimensions 1 and 2, indexed by dimension 3.
#' @param row.mean By default, \code{FALSE}. If \code{TRUE}, will fit a
#'    common mean within each row. If both this and \code{col.mean} are
#'    \code{TRUE}, there will be a common mean for the entire matrix.
#' @param col.mean By default, \code{FALSE}. If \code{TRUE}, will fit a
#'    common mean within each row. If both this and \code{row.mean} are
#'    \code{TRUE}, there will be a common mean for the entire matrix.
#' @param row.variance Imposes a variance structure on the rows. Either
#'    'none', 'AR(1)', 'CS' for 'compound symmetry', 'Correlation' for a
#'    correlation matrix, or 'Independence' for
#'    independent and identical variance across the rows.
#'    Only positive correlations are allowed for AR(1) and CS covariances.
#'    Note that while maximum likelihood estimators are available (and used) for
#'    the unconstrained variance matrices, \code{optim} is used for any
#'    constraints so it may be considerably slower.
#' @param col.variance  Imposes a variance structure on the columns.
#'     Either 'none', 'AR(1)', 'CS', 'Correlation', or 'Independence'.
#'     Only positive correlations are allowed for
#'     AR(1) and CS.
#' @param tol Convergence criterion. Measured against square deviation
#'    between iterations of the two variance-covariance matrices.
#' @param max.iter Maximum possible iterations of the algorithm.
#' @param U (optional) Can provide a starting point for the \code{U} matrix.
#'    By default, an identity matrix.
#' @param V (optional) Can provide a starting point for the \code{V} matrix.
#'    By default, an identity matrix.
#' @param ... (optional) additional arguments can be passed to \code{optim}
#'    if using restrictions on the variance.
#'
#' @return Returns a list with a the following elements:
#' \describe{
#'     \item{\code{mean}}{the mean matrix}
#'     \item{\code{scaling}}{the scalar variance parameter
#'          (the first entry of the covariances are restricted to unity)}
#'     \item{\code{U}}{the between-row covariance matrix}
#'     \item{\code{V}}{the between-column covariance matrix}
#'     \item{\code{iter}}{the number of iterations}
#'     \item{\code{tol}}{the squared difference between iterations of
#'          the variance matrices at the time of stopping}
#'     \item{\code{logLik}}{vector of log likelihoods at each iteration.}
#'     \item{\code{convergence}}{a convergence flag, \code{TRUE} if converged.}
#'     \item{\code{call}}{The (matched) function call.}
#'    }
#'
#' @references
#' Pierre Dutilleul.  The MLE algorithm for the matrix normal
#' distribution.
#'     Journal of Statistical Computation and Simulation, (64):105–123, 1999.
#'
#'     Gupta, Arjun K, and Daya K Nagar. 1999. Matrix Variate Distributions.
#'     Vol. 104. CRC Press. ISBN:978-1584880462
#' @export
#' @seealso \code{\link{rmatrixnorm}} and \code{\link{MLmatrixt}}
#'
#' @examples
#' set.seed(20180202)
#' # simulating from a given density
#' A <- rmatrixnorm(
#'   n = 100, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
#'   L = matrix(c(2, 1, 0, .1), nrow = 2), list = TRUE
#' )
#' # finding the parameters by ML estimation
#' results <- MLmatrixnorm(A, tol = 1e-5)
#' print(results)
MLmatrixnorm <- function(data, row.mean = FALSE, col.mean = FALSE,
                         row.variance = "none", col.variance = "none",
                         tol = 10 * .Machine$double.eps^0.5, max.iter = 100,
                         U, V, ...) {
  if (class(data) == "list") {
    data <- array(unlist(data),
      dim = c(
        nrow(data[[1]]),
        ncol(data[[1]]), length(data)
      )
    )
  }
  if (!all(
    is.numeric(data), is.numeric(tol),
    is.numeric(max.iter)
  )) {
    stop("Non-numeric input. ")
  }
  if (!(missing(U))) {
    if (!(is.numeric(U))) stop("Non-numeric input.")
  }
  if (!(missing(V))) {
    if (!(is.numeric(V))) stop("Non-numeric input.")
  }
  row_set_var <- FALSE


  rowvarparse <- .varparse(row.variance)
  row_set_var <- rowvarparse$varflag
  row.variance <- rowvarparse$varopt
  col_set_var <- FALSE



  colvarparse <- .varparse(col.variance)
  col_set_var <- colvarparse$varflag
  col.variance <- colvarparse$varopt

  dims <- dim(data)

  if (max(dims[1] / dims[2], dims[2] / dims[1]) > (dims[3] - 1)) {
    stop("Need more observations to estimate parameters.")
  }
  # don't have initial starting point for U and V, start with diag.
  if (missing(U)) {
    U <- diag(dims[1])
  }
  if (missing(V)) {
    V <- diag(dims[2])
  }
  # mu <- apply(data, c(1, 2), mean)
  mu <- rowMeans(data, dims = 2)
  if (row.mean) {
    # make it so that the mean is constant within a row
    mu <- matrix(rowMeans(mu), nrow = dims[1], ncol = dims[2])
  }
  if (col.mean) {
    # make it so that the mean is constant within a column
    mu <- matrix(colMeans(mu), nrow = dims[1], ncol = dims[2], byrow = TRUE)
  }
  # if both are true, this makes it so the mean is constant all over
  swept_data <- sweep(data, c(1, 2), mu)
  iter <- 0
  error.term <- 1e+40
  if (col_set_var) {
    if (V[1, 2] > 0) {
      rho_col <- V[1, 2] / max(V[1, 1], 1)
    } else {
      inter_v <- txax(swept_data, .5 * (U + t(U)))
      V <- rowSums(inter_v, dims = 2) / (dims[3] * dims[1])
      if (col.variance == "AR(1)") {
        V <- stats::cov2cor(V)
        rho_col <- V[1, 2] / V[1, 1]
      }
      if (col.variance == "CS") {
        V <- stats::cov2cor(V)
        rho_col <- mean(V[1, ] / V[1, 1])
      }
      if (col.variance == "I") rho_col <- 0
      if (rho_col > .9) rho_col <- .9
      if (rho_col < 0) rho_col <- 0
      V <- varmatgenerate(dims[2], rho_col, col.variance)
    }
  }

  if (row_set_var) {
    if (U[1, 2] > 0) {
      rho_row <- U[1, 2] / max(U[1, 1], 1)
    } else {
      inter_u <- xatx(swept_data, 0.5 * (V + t(V)))
      U <- rowSums(inter_u, dims = 2) / (dims[3] * dims[2])
      if (row.variance == "AR(1)") {
        U <- stats::cov2cor(U)
        rho_row <- U[1, 2] / U[1, 1]
      }
      if (row.variance == "CS") {
        U <- stats::cov2cor(U)
        rho_row <- mean(U[1, ] / U[1, 1])
      }
      if (row.variance == "I") rho_row <- 0
      if (rho_row > .9) rho_row <- .9
      if (rho_row < 0) rho_row <- 0
      # if (row.variance == "cor") U = stats::cov2cor(U) else
      U <- varmatgenerate(dims[1], rho_row, row.variance)
    }
  }

  varflag <- FALSE
  log_lik_vec <- numeric(0)
  while (iter < max.iter && error.term > tol && (!varflag)) {

    # make intermediate matrix, then collapse to final version
    if (col_set_var) {
      var <- V[1, 1]
      invV <- chol2inv(chol.default(V / var))
      invU <- chol2inv(chol.default((U)))
      var <- sum(apply(
        matrix(swept_data, ncol = dims[3]), 2,
        # function(x) txax(x, invV %x% invU))) / (prod(dims))
        # not sure why txax() doesn't work here
        function(x) {
          crossprod(
            x,
            invV %x% invU
          ) %*% x
        }
      )) / (prod(dims))
      if (col.variance != "I") {
        tmp <- txax(swept_data, 0.5 * (U + t(U)))
        tmpsummary <- matrix(rowSums(tmp, FALSE, dims = 2), nrow = dims[2])
        nLL <- function(theta) {
          Vmat <- varinv(dims[2], theta, TRUE, col.variance) / var # try it
          B <- Vmat %*% tmpsummary
          # solved derivative, need to find where this is zero:
          0.5 * dims[1] * dims[3] * vardet(dims[2], theta, TRUE, col.variance) -
            (.5) * sum(diag(B)) # problem was wrong constant
        }
        if (!isTRUE(sign(nLL(0)) * sign(nLL(.999)) <= 0)) {
          warning("Endpoints of derivative of likelihood do not have
                   opposite sign. Check variance specification.")
          rho_col <- 0
          varflag <- TRUE
        } else {
          fit0 <- stats::uniroot(nLL, c(0, .999), ...)
          rho_col <- fit0$root
        }
      }
      new_v <- var * varmatgenerate(dims[2], rho_col, col.variance)
    } else {
      inter_v <- txax(swept_data, 0.5 * (U + t(U)))
      new_v <- rowSums(inter_v, dims = 2) / (dims[3] * dims[1])
      if (col.variance == "cor") {
        vartmp <- exp(mean(log(diag(new_v)))) # matrix should be pos definite
        if (!is.finite(vartmp)) {
          vartmp <- 1
          varflag <- TRUE
          warning("Variance estimate for correlation matrix
                        not positive definite.")
        }
        new_v <- vartmp * stats::cov2cor(new_v)
      }
    }
    if (row.variance == "I") {
      new_u <- diag(dims[1])
    } else if (row_set_var) {
      tmp <- xatx(swept_data, 0.5 * (V + t(V)))
      tmpsummary <- matrix(rowSums(tmp, dims = 2), nrow = dims[1])
      nLL <- function(theta) {
        Umat <- varinv(dims[1], theta, TRUE, row.variance)
        B <- Umat %*% tmpsummary
        # solved derivative, need to find where this is zero:
        0.5 * dims[2] * dims[3] * vardet(dims[1], theta, TRUE, row.variance) -
          (.5) * sum(diag(B)) # problem was wrong constant
      }
      if (!isTRUE(sign(nLL(0)) * sign(nLL(.999)) <= 0)) {
        warning("Endpoints of derivative of likelihood do not have
                 opposite sign.  Check variance specification.")
        rho_row <- 0
        varflag <- TRUE
      } else {
        fit0 <- stats::uniroot(nLL, c(0, .999), ...)
        rho_row <- fit0$root
      }
      new_u <- varmatgenerate(dims[1], rho_row, row.variance)
    } else {
      inter_u <- xatx(swept_data, 0.5 * (new_v + t(new_v)))
      new_u <- rowSums(inter_u, dims = 2) / (dims[3] * dims[2])
      new_u <- new_u / (new_u[1, 1])
      if (row.variance == "cor") new_u <- stats::cov2cor(new_u)
    }
    # only identifiable up to a constant, so have to fix something at 1

    error.term <- sum((new_v - V)^2) + sum((new_u - U)^2)
    V <- new_v
    U <- new_u
    log_lik <- sum(dmatrixnorm(data, mu, U = U, V = V, log = TRUE))
    log_lik_vec <- c(log_lik_vec, log_lik)
    iter <- iter + 1
  }
  if (iter >= max.iter || error.term > tol || varflag) {
    warning("Failed to converge")
  }

  converged <- !(iter >= max.iter || error.term > tol || varflag)
  log_lik <- 0

  log_lik <- sum(dmatrixnorm(data, mu, U = U, V = V, log = TRUE))
  return(list(
    mean = mu,
    U = U,
    V = V / V[1, 1],
    var = V[1, 1],
    iter = iter,
    tol = error.term,
    logLik = log_lik_vec,
    convergence = converged,
    call = match.call()
  ))
}
