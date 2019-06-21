#' Distribution functions for the matrix variate t distribution.
#' @family matrix-variate
#'
#' Density and random generation for the matrix variate t distribution.
#'
#' @param n number of observations for generation
#' @param x quantile for density
#' @param df  degrees of freedom (\eqn{>0}, may be non-integer),
#'    \code{df = 0, Inf} is allowed and will return a normal distribution.
#' @param mean \eqn{p \times q}{p * q} This is really a 'shift' rather than a mean,
#'    though the expected value will be equal to this if
#'    \eqn{df > 2}
#' @param L \eqn{p \times p}{p * p}  matrix specifying relations among the rows.
#'     By default, an identity matrix.
#' @param R \eqn{q \times q}{q * q}  matrix specifying relations among the
#'     columns. By default, an identity matrix.
#' @param U \eqn{LL^T}  - \eqn{p \times p}{p * p}   positive definite matrix for
#'    rows, computed from \eqn{L} if not specified.
#' @param V \eqn{R^T R}  - \eqn{q \times q}{q * q}  positive definite matrix for
#'    columns, computed from \eqn{R}  if not specified.
#' @param list Defaults to \code{FALSE} . If this is \code{TRUE} , then the
#'    output will be a list of matrices.
#' @param array If \eqn{n = 1}  and this is not specified and \code{list}  is
#'    \code{FALSE} , the function will return a matrix containing the one
#'    observation. If \eqn{n > 1} , should be the opposite of \code{list} .
#'    If \code{list}  is \code{TRUE} , this will be ignored.
#' @param force In \code{rmatrix}: if \code{TRUE}, will take the input of
#'    \code{R} directly - otherwise uses \code{V} and uses Cholesky
#'    decompositions. Useful for generating degenerate t-distributions.
#'    Will also override concerns about potentially singular matrices
#'    unless they are not, in fact, invertible.
#' @param log logical; in \code{dmatrixt}, if \code{TRUE}, probabilities
#'    \code{p} are given as \code{log(p)}.
#' @return \code{rmatrixt} returns either a list of \eqn{n}
#'    \eqn{p \times q}{p * q}  matrices ora \eqn{p \times q \times n}{p * q * n}
#'    array.
#'
#'    \code{dmatrixt} returns the density at quantile \code{x}.
#' @details
#' The matrix \eqn{t}-distribution is parameterized slightly
#'  differently from the univariate and multivariate \eqn{t}-distributions
#'  - the variance is scaled by a factor of \code{1/df}.
#'  In this parameterization, the variance for a \eqn{1 \times 1}{1 * 1} matrix variate
#'  \eqn{t}-distributed random variable with identity variance matrices is
#'  \eqn{1/(df-2)} instead of \eqn{df/(df-2)}. A Central Limit Theorem
#'  for the matrix variate \eqn{T} is then that as \code{df} goes to
#'  infinity, \eqn{MVT(0, df, I_p, df*I_q)} converges to
#'  \eqn{MVN(0,I_p,I_q)}.
#'
#' @seealso \code{\link{rt}} and \code{\link[stats]{Distributions}}
#'
#' @export
#'
#' @examples
#' set.seed(20180202)
#' rmatrixt(n=1,df=10,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'    L=matrix(c(2,1,0,.1),nrow=2),list=FALSE)
#' summary(rt(n=100,df=10))
#' summary(rmatrixt(n=100,df=10,matrix(0)))
#' set.seed(20180204)
#' x = rmatrixt(n=1,mean=matrix(0),df=1)
#' dt(x,1)
#' dmatrixt(x,df=1)
#'

rmatrixt <- function(n, df, mean,
                     L = diag(dim(as.matrix(mean))[1]),
                     R = diag(dim(as.matrix(mean))[2]),
                     U = L %*% t(L),
                     V = t(R) %*% R,
                     list = FALSE,
                     array = NULL,
                     force = FALSE) {
  if (!(all(is.numeric(n),is.numeric(df), is.numeric(mean),
            is.numeric(L), is.numeric(R),
            is.numeric(U),is.numeric(V)))) stop("Non-numeric input. ")
  if (!(n > 0))
    stop("n must be > 0. n =", n)

  if (length(df) != 1) stop("Length of df must be 1. length = ", length(df))

  if (((is.null(df)) || is.na(df) || (df < 0)))
    stop("df must be >= 0. df =", df)
  if (df == 0 || is.infinite(df))
    return(rmatrixnorm(n = n, mean = mean, U = U,
                       V = V, list = list, array = array))
  mean <- as.matrix(mean)
  U <- as.matrix(U)
  V <- as.matrix(V)
  if (!symm.check(U)) stop("U not symmetric.")
  if (!symm.check(V)) stop("V not symmetric.")
  dims <- dim(mean)

  if (!(dims[1] == dim(U)[2] && dim(U)[1] == dim(U)[2] &&
        dims[2] == dim(V)[1] && dim(V)[1] == dim(V)[2])) {
    stop("Non-conforming dimensions.", dims, dim(U), dim(V))
  }
  if (force && !missing(R)) cholV <- R else cholV <- chol(V)

  if (min(diag(cholV)) < 1e-6 && !force) {
    stop("Potentially singular covariance, use force = TRUE if intended. ",
         min(diag(cholV)))
  }

  nobs <- prod(dims)*n
  mat <- array(stats::rnorm(nobs), dim = c(dims,n))

  cholU <- CholWishart::rInvCholWishart(n, df + dims[1] - 1,U)

  result <- array(dim = c(dims,n))

  for (i in seq(n)) {
    result[ , , i] <- mean + (crossprod(cholU[ , , i], mat[ , , i])) %*% (cholV)
  }

  if (n == 1 && list == FALSE && is.null(array)) {
    return(array(result[,,1], dim = dims[1:2]))
    # if n = 1 and you don't specify arguments, it just returns a matrix
  }
  if (list) {
    return(lapply(seq(dim(result)[3]), function(x) result[ , , x]))
  }
  if (!(list) && !(is.null(array))) {
    if (!(array))
      warning("list FALSE and array FALSE, returning array")
  }
  return(result)
}

#' @describeIn rmatrixt Density function for matrix variate t distribution
#' @export
dmatrixt <- function(x, df, mean = matrix(0, p, n),
           L = diag(p),
           R = diag(n), U = L %*% t(L),
           V = t(R) %*% R, log = FALSE){

  dims <- dim(x)
  if (is.null(dims) || length(dims) == 1) x <- matrix(x)

  dims <- dim(x)
  if (length(dims) == 2) x <- array(x, dim = (dims <- c(dims,1)))
  p <- dims[1]
  n <- dims[2]
  if (!(all(is.numeric(x), is.numeric(mean), is.numeric(L), is.numeric(R),
            is.numeric(U),is.numeric(V)))) stop("Non-numeric input. ")
  if (length(df) != 1) stop("Length of df must be 1. length = ", length(df))
  if ( ((is.null(df)) || is.na(df) || (df < 0)))
    stop("df must be >= 0. df =", df)

  mean <- as.matrix(mean)
  U <- as.matrix(U)
  V <- as.matrix(V)
  if (!symm.check(U)) stop("U not symmetric.")
  if (!symm.check(V)) stop("V not symmetric.")

  if (!(dims[1] == dim(U)[2] && dim(U)[1] == dim(U)[2] &&
        dims[2] == dim(V)[1] && dim(V)[1] == dim(V)[2])) {
    stop("Non-conforming dimensions.", dims, dim(U),dim(V))
  }
  if ( (df == 0 || is.infinite(df)) )
    return(dmatrixnorm(x, mean = mean, U = U, V = V, log = log))

  # gammas is constant
  # this could be shifted into C++ but I don't want to pull out of CholWishart
  gammas <- as.numeric(CholWishart::lmvgamma((0.5) * (df + dims[1] + dims[2] - 1), dims[1]) -
    0.5 * dims[1]*dims[2] * log(pi) -
    CholWishart::lmvgamma(0.5 * (df + dims[1] - 1), dims[1]))

  results = as.numeric(dmat_t_calc(x, df, mean, U, V))
  results = results + gammas
  if (log) {
    return(results)
  } else {
    return(exp(results))
  }
}


#' Distribution functions for matrix variate inverted t distributions
#'
#' Generate random draws from the inverted matrix
#'    variate t distribution
#' @family matrix-variate
#' @inheritParams rmatrixt
#' @return \code{rmatrixinvt} returns either a list of \eqn{n}
#'    \eqn{p \times q}{p * q}  matrices or
#'    a \eqn{p \times q \times n}{p * q * n}  array.
#'
#'    \code{dmatrixinvt} returns the density at quantile \code{x}.

#' @export
#' @examples
#' A<-rmatrixinvt(n = 2, df = 10, diag(4))
#' dmatrixinvt(A[,,1], df = 10, mean = diag(4))
rmatrixinvt <- function(n, df, mean,
                        L = diag(dim(as.matrix(mean))[1]),
                        R = diag(dim(as.matrix(mean))[2]),
                        U = L %*% t(L), V = t(R) %*% R,
                        list=FALSE, array = NULL) {
  if (!(all(is.numeric(df), is.numeric(mean), is.numeric(L), is.numeric(R),
            is.numeric(U),is.numeric(V)))) stop("Non-numeric input. ")
  if (length(df) != 1) stop("Length of df must be 1. length = ", length(df))
  if ( ((is.null(df)) || is.na(df) || (df < 0)))
    stop("df must be >= 0. df =", df)
  if (!(n > 0))
    stop("n must be > 0.", n)
  mean <- as.matrix(mean)
  U <- as.matrix(U)
  V <- as.matrix(V)
  if (!symm.check(U)) stop("U not symmetric.")
  if (!symm.check(V)) stop("V not symmetric.")
  dims <- dim(mean)
  # checks for conformable matrix dimensions
  if (!(dims[1] == dim(U)[2] && dim(U)[1] == dim(U)[2] &&
        dims[2] == dim(V)[1] && dim(V)[1] == dim(V)[2])) {
    stop("Non-conforming dimensions.", dims, dim(U), dim(V))
  }

  nobs <- prod(dims)*n
  mat <- array(stats::rnorm(nobs), dim = c(dims,n))

   S <- stats::rWishart(n, df + dims[1] - 1, diag(dims[1]))

  result = rmat_inv_t_calc(S, mat, U, V, mean)


  if (n == 1 && list == FALSE && is.null(array)) {
    return(array(result[,,1], dim = dims[1:2]))
    # if n = 1 and you don't specify arguments, it just returns a matrix
  }
  if (list) {
    return(lapply(seq(dim(result)[3]), function(x) result[ , , x]))
  }
  if (!(list) && !(is.null(array))) {
    if (!(array))
      warning("list FALSE and array FALSE, returning array")
  }
  return(result)
}



#' @describeIn rmatrixinvt Density function for inverted matrix t-distributions.
#' @export
dmatrixinvt <- function(x, df, mean = matrix(0, p, n),
                     L = diag(p),
                     R = diag(n), U = L %*% t(L),
                     V = t(R) %*% R, log = FALSE){

  dims <- dim(x)
  if (is.null(dims) || length(dims) == 1) x <- matrix(x)

  dims <- dim(x)
  if (length(dims) == 2) x <- array(x, dim = (dims <- c(dims,1)))
  p <- dims[1]
  n <- dims[2]
  if (!(all(is.numeric(x), is.numeric(mean), is.numeric(L), is.numeric(R),
            is.numeric(U),is.numeric(V)))) stop("Non-numeric input. ")
  if (length(df) != 1) stop("Length of df must be 1. length = ", length(df))
  if ( ((is.null(df)) || is.na(df) || (df < 0)))
    stop("df must be >= 0. df =", df)

  mean <- as.matrix(mean)
  U <- as.matrix(U)
  V <- as.matrix(V)
  if (!symm.check(U)) stop("U not symmetric.")
  if (!symm.check(V)) stop("V not symmetric.")

  if (!(dims[1] == dim(U)[2] && dim(U)[1] == dim(U)[2] &&
        dims[2] == dim(V)[1] && dim(V)[1] == dim(V)[2])) {
    stop("Non-conforming dimensions.", dims, dim(U),dim(V))
  }
  gammas <- as.numeric(CholWishart::lmvgamma((0.5) * (df + dims[1] + dims[2] - 1), dims[1]) -
    0.5 * prod(dims[1:2]) * log(pi) - CholWishart::lmvgamma(0.5 * (df + dims[1] - 1), dims[1]))

    results = as.numeric(dmat_inv_t_calc(x, df,  mean, U, V))
    if( any(is.nan(results))) warning("warning: probability distribution undefined when det < 0.")
  results <- gammas + results
  if (log) {
    return(results)
  } else {
    return(exp(results))
  }
}
