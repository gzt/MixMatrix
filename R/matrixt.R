#' Distribution functions for the matrix variate t-distribution.
#' @family matrixt
#'
#' Density and random generation for the matrix variate t-distribution.
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
#' @param R \eqn{q \times q}{q * q}  matrix specifying relations among the columns.
#'     By default, an identity matrix.
#' @param U \eqn{LL^T}  - \eqn{p \times p}{p * p}   positive definite matrix for rows,
#'    computed from \eqn{L} if not specified.
#' @param V \eqn{R^T R}  - \eqn{q \times q}{q * q}  positive definite matrix for columns,
#'    computed from \eqn{R}  if not specified.
#' @param list Defaults to \code{FALSE} . If this is \code{TRUE} , then the
#'    output will be a list of matrices.
#' @param array If \eqn{n = 1}  and this is not specified and \code{list}  is
#'    \code{FALSE} , the function will return a matrix containing the one
#'    observation. If \eqn{n > 1} , should be the opposite of \code{list} .
#'    If \code{list}  is \code{TRUE} , this will be ignored.
#' @param force In \code{rmatrix}: if TRUE, will take the input of \code{R}
#'    directly - otherwise uses \code{V} and uses Cholesky
#'    decompositions. Useful for generating degenerate t-distributions.
#'    Will also override concerns about potentially singular matrices
#'    unless they are not, in fact, invertible.
#' @param log logical; in \code{dmatrixt}, if TRUE, probabilities p are given as log(p).
#' @return \code{rmatrixt} returns either a list of \eqn{n}  \eqn{p \times q}{p * q}  matrices or
#'    a \eqn{p \times q \times n}{p * q * n}  array.
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
  # should probably do better error checking, checks for
  # conformable matrix dimensions
  if (!(dims[1] == dim(U)[2] && dim(U)[1] == dim(U)[2] &&
        dims[2] == dim(V)[1] && dim(V)[1] == dim(V)[2])) {
    stop("Non-conforming dimensions.", dims, dim(U), dim(V))
  }
  if (force && !missing(R)) cholV <- R else cholV <- chol(V)

  if (min(diag(cholV)) < 1e-6 && !force) {
    stop("Potentially singular covariance, use force = TRUE if intended. ",
         min(diag(cholV)))
  }
  # this is not necessary, the correct input is
  # U, not U^-1
  #solveU = chol2inv(chol.default(U))

  nobs <- prod(dims)*n
  mat <- array(stats::rnorm(nobs), dim = c(dims,n))

  # USigma <- stats::rWishart(n, df + dims[1] - 1, (1/df) * solveU)

  cholU <- CholWishart::rInvCholWishart(n, df + dims[1] - 1,U)

  result <- array(dim = c(dims,n))
  for (i in seq(n)) {
    # just changed: pretty sure this is supposed to be t(chol(...))
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

#' @describeIn rmatrixt Density function for matrix variate t-distribution
#' @export
dmatrixt <- function(x, df, mean = matrix(0, p, n),
           L = diag(p),
           R = diag(n), U = L %*% t(L),
           V = t(R) %*% R, log = FALSE){
  # x <- as.matrix(x)
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
  gammas <- CholWishart::lmvgamma((0.5) * (df + dims[1] + dims[2] - 1), dims[1]) -
    0.5 * dims[1]*dims[2] * log(pi) -
    CholWishart::lmvgamma(0.5 * (df + dims[1] - 1), dims[1])

  ## starting here, can try to put in CPP

  xm <- array(0, dim = c(p,n,dims[3]))
  for(i in seq(dims[3])){
  xm[, , i] <- x[, , i] - mean
  }
  cholU <- chol.default(U)
  cholV <- chol.default(V)
  # breaking equation into two parts: the integrating constants (gammas)
  # and the matrix algebra parts (mats) done on the log scale

  if (any(diag(cholU) < 1e-6) || any(diag(cholV) < 1e-6)) stop("non-invertible matrix", min(diag(cholU)), min(diag(cholV)))
  logdetU <- 2*sum(log(diag(cholU)))
  logdetV <- 2*sum(log(diag(cholV)))
  results <- rep(gammas, dims[3])
  for(i in seq(dims[3])){
  m <- diag(dims[1]) + crossprod(chol2inv(cholU), matrix(xm[,,i], dims[1], dims[2])) %*%
    tcrossprod(chol2inv(cholV),matrix(xm[,,i], dims[1], dims[2]))
  #m <- diag(dims[1]) + chol2inv(cholU) %*% xatx(xm, chol2inv(cholV))
  # - 0.5 * dims[2] * dims[1]*log(df) term disappears
  mats <- -0.5 * dims[2] * (logdetU)  -
    0.5 * dims[1] * logdetV -
    0.5 * (df + dims[1] + dims[2] - 1) * log(det(m))

  results[i] <- results[i] + mats
  }
  if (log) {
    return(results)
  } else {
    return(exp(results))
  }
}


#' Distribution functions for matrix variate inverted t-distributions
#'
#' Generate random draws from the inverted matrix
#'    variate t-distribution
#' @family matrixt
#' @inheritParams rmatrixt
#' @return \code{rmatrixinvt} returns either a list of \eqn{n}  \eqn{p \times q}{p * q}  matrices or
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

  Usqrt <- posmatsqrt(U)
  Vsqrt <- posmatsqrt(V)
  if (any(!is.finite(Usqrt)) || any(!is.finite(Vsqrt))) stop("Non-positive definite U or V matrix.")
  nobs <- prod(dims)*n
  mat <- array(stats::rnorm(nobs), dim = c(dims,n))

  S <- stats::rWishart(n, df + dims[1] - 1, diag(dims[1]))

  SXX <- array(0,dim = c(dims[1],dims[1],n))

  for (i in seq(n)) {
    # if there's a way to do this with apply I want to see it
    SXX[ , , i] <- array((posmatsqrtinv(S[ , , i] +
                                          tcrossprod(mat[ , , i]))),
                         dim = c(dims[1],dims[1]))
  }

  result <- array(dim = c(dims,n))

  for (i in seq(n)) {
    # if there's a way to do this with apply I want to see it
    result[ , , i] <- Usqrt %*% SXX[ , , i] %*% mat[ , , i] %*% Vsqrt + mean
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



#' @describeIn rmatrixinvt Density function for inverted matrix t-distributions.
#' @export
dmatrixinvt <- function(x, df, mean = array(0, dim(as.matrix(x))[1:2]),
                        L = diag(dim(as.matrix(x))[1]),
                        R = diag(dim(as.matrix(x))[2]),
                        U = L %*% t(L), V = t(R) %*% R, log = FALSE) {
  if (!(all(is.numeric(x),is.numeric(df), is.numeric(mean), is.numeric(L),
            is.numeric(R), is.numeric(U),
            is.numeric(V)))) stop("Non-numeric input. ")
  if (length(df) != 1) stop("Length of df must be 1. length = ", length(df))
  if ( ((is.null(df)) || is.na(df) || (df < 0)))
    stop("df must be >= 0. df =", df)
  if (df <= 0 || is.infinite(df)) stop("Invalid input for df,
                                        must have 0 < df < Inf, df = ", df)
  x <- as.matrix(x)
  mean <- as.matrix(mean)
  U <- as.matrix(U)
  V <- as.matrix(V)
  if (!symm.check(U)) stop("U not symmetric.")
  if (!symm.check(V)) stop("V not symmetric.")
  dims <- dim(x)
  if (!(dims[1] == dim(U)[2] && dim(U)[1] == dim(U)[2] &&
        dims[2] == dim(V)[1] && dim(V)[1] == dim(V)[2])) {
    stop("Non-conforming dimensions.", dims, dim(U), dim(V))
  }

  gammas <- CholWishart::lmvgamma((0.5) * (df + sum(dims) - 1), dims[1]) -
    0.5 * prod(dims) * log(pi) - CholWishart::lmvgamma(0.5 * (df + dims[1] - 1), dims[1])

  ## starting here can put in CPP
  xm <- x - mean
  cholU <- chol.default(U)
  cholV <- chol.default(V)

  if ( (any(diag(cholU) < 1e-6) || any(diag(cholV) < 1e-6))) stop("non-invertible matrix", min(diag(cholU)), min(diag(cholV)))
  logdetU <- 2*sum(log(diag(cholU)))
  logdetV <- 2*sum(log(diag(cholV)))
  # breaking equation into two parts: the integrating constants
  # (gammas) and the matrix algebra parts (mats) done on the log scale
  # note I did not do the DF correction as for the matrix t distribution



  matrixterms <- diag(dims[1]) -
    #chol2inv(cholU) %*% xatx(xm,  chol2inv(cholV))
    chol2inv(cholU) %*% xm %*% tcrossprod(chol2inv(cholV),(xm))

  mats <- -0.5 * dims[2] * logdetU - 0.5 * dims[1] * logdetV -
    0.5 * (df - 2) * log(det(matrixterms))
  results <- gammas + mats
  if (log) {
    return(results)
  } else {
    return(exp(results))
  }
}
