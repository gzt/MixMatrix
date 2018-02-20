
#' Generate a unit AR(1) covariance matrix
#'
#' @description generate AR(1) correlation matrices
#'
#' @param n number of columns/rows
#' @param rho correlation parameter
#'
#' @return Toeplitz \eqn{n x n} matrix with 1 on the diagonal and \eqn{rho^k} on
#'    the other diagonals, where \eqn{k} is distance from the main diagonal.
#'    Used internally but it is useful for generating your own random matrices.
#' @seealso \code{\link[stats]{toeplitz}}
#' @export
#'
#' @examples
#' toepgenerate(6,.9)
toepgenerate <- function(n, rho) {
  if (n <= 1)
    stop("n must be greater than 1.")
  if (rho >= 1)
    stop("rho must be a correlation less than 1.")
  if (rho <= -1)
    stop("rho must be a correlation greater than -1.")
  if (rho < 0)
    warning("Rho = ", rho, " and should be greater than 0.")
  if (rho > 0.99)
    warning("Rho = ", rho, " high correlation may cause numerical problems.")
  X <- stats::toeplitz(c(1, rho^(1:(n - 1))))
  return(X)
}


#' Generate a compound symmetric correlation matrix
#'
#' @param n number of dimensions
#' @param rho off-diagonal element - a correlation between -1 and 1. Will warn if less than 0.
#'
#' @return returns an \eqn{n x n} matrix with 1 on the diagonal and \code{rho} on the off-diagonal.
#' @export
#'
#' @examples
#' CSgenerate(3,.5)

CSgenerate <- function(n,rho) {
  if (n <= 1)
    stop("n must be greater than 1.")
  if (rho >= 1)
    stop("rho must be a correlation less than 1.")
  if (rho <= -1)
    stop("rho must be a correlation greater than -1.")
  if (rho < 0)
    warning("Rho = ", rho, " and should be greater than 0.")
  if (rho > 0.99)
    warning("Rho = ", rho, " high correlation may cause numerical problems.")
  A <- matrix(rho, nrow = n,ncol = n)
  diag(A) <- 1
  A
}


#' Quick symmetry check
#'
#' Quick check whether matrix input is symmetric -
#' checks sum of absolute differences of transposes
#' as well as dimensions. Not robust, so only an
#' internal function to be used with known safe input.
#'
#' @param A Numeric real matrix. Does not check if real.
#' @param tol tolerance - note that if you have a big matrix
#'    it may need to be specified as it's a sum of entries.
#'
#' @return logical TRUE if symmetric FALSE otherwise.
#' Not as robust as \code{isSymmetric()}.
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' A = toepgenerate(5,.9)
#' symm.check(A)
#' A[1,2] = 5
#' symm.check(A)}
symm.check <- function(A, tol = 10 * (.Machine$double.eps)^.5) {
  # if (!is.matrix(A)) return(FALSE)
  # if (!is.numeric(A)) return(FALSE)
  # commented those out because it is always checked before running symm.check.
  dims <- dim(A)
  if (dims[1] != dims[2]) {
    return(FALSE)
  }
  return(sum(abs(A - t(A))) < prod(dims)*tol)
}


#' Select a variance structure to generate.
#'
#' @param n number of dimensions
#' @param rho parameter for selected variance structure.
#' @param variance variance structure - AR(1) or CS.
#'
#' @return Specified matrix structure
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{varmatgenerate(5,.2,"AR(1)")}
#'
varmatgenerate <- function(n, rho, variance) {
  if (variance == "I" || variance == "independence" || variance == "Independence") variance = "I"
  if (variance == "AR(1)") return(toepgenerate(n,rho))
  if (variance == "CS") return(CSgenerate(n,rho))
  if (variance == "I") return(diag(n))
  else stop("Bad covariance structure input.", variance)

}

#' Determinant selector for chosen covariance matrix.
#'
#' @param n dimensions
#' @param rho off-diagonal parameter
#' @param deriv logical whether to return the determinant or the derivative of
#'     the log of the determinant
#' @param variance  variance structure - AR(1) or CS.
#'
#' @return Determinant or derivative of log-inverse for the specified matrix structure.
#' @keywords internal
#'
#' @examples
#' \dontrun{vardet(5,.5,TRUE, "AR(1)"}
vardet <- function(n, rho, deriv, variance){

  if (variance == "AR(1)") return(ARdet(n,rho, deriv))

  if (variance == "CS") return(CSdet(n,rho, deriv))

  else stop("Bad covariance structure input.", variance)
}

#' Inverse selector for chosen covariance matrix.
#'
#' @param n dimensions
#' @param rho off-diagonal parameter
#' @param deriv logical whether to return the inverse or the derivative of
#'     the inverse
#' @param variance  variance structure - AR(1) or CS.
#'
#' @return The inverse or derivative of the inverse of the selected matrix structure.
#' @keywords internal
#'
#' @examples
#' \dontrun{varinv(5,.5,TRUE, "AR(1)"}
varinv <- function(n, rho, deriv, variance){

  if (variance == "AR(1)") return(invAR(n,rho, deriv))

  if (variance == "CS") return(invCS(n,rho, deriv))

  else stop("Bad covariance structure input. ", variance)
}


#' Determinant for an AR(1) covariance matrix.
#'
#' @param n dimensions
#' @param rho off-diagonal parameter
#' @param deriv logical whether to return the determinant or the derivative of
#'     the log of the determinant
#'
#' @return determinant of an AR(1) covariance matrix. If \code{deriv} is specified,
#'     will return the derivative of \eqn{\log |\Sigma^{-1}|}.
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' ARdet(toepgenerate(5,.4,deriv=TRUE))
#' }
#'
ARdet <- function(n, rho, deriv = FALSE) {

  if (!deriv) return((1 - rho^2)^(n - 1))
  else (1 - n) * (-2 * rho) / (1 - rho^2)

}


#' Determinant for an CS covariance matrix.
#'
#' @param n dimensions
#' @param rho off-diagonal parameter
#' @param deriv logical whether to return the determinant or the
#' derivative of the log determinant of the inverse
#'
#' @return determinant of an CS covariance matrix
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' CSdet(toepgenerate(5,.4,deriv=TRUE))
#' }
#'
CSdet <- function(n, rho, deriv = FALSE) {

  if (!deriv) return((1 + rho*(n - 1))*(1 - rho)^(n - 1))
  else (n - 1) * n * rho /((1 - rho) * ((n - 1) * rho + 1) )
  # -n*(n - 1) * rho * (1 - rho)^(n - 2)
}


#' Returns the inverse of an AR(1) covariance matrix or its derivative
#'
#' @param n dimensions of matrix
#' @param rho correlation parameter
#' @param deriv logical. if TRUE will output the derivative of the inverse matrix.
#'
#' @return Matrix of the inverse of the AR(1) covariance matrix (or its inverse)
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' invAR(5,.5, TRUE)
#' }
invAR <- function(n,rho, deriv = FALSE){
  if (!(n > 1)) stop("n needs to be greater than 1")
  if (!(rho < 1 && rho > -1)) stop("rho needs to be < 1")
  if (deriv) {
    #-(rho^2+1)/(1-rho^2)^2  , 4*rho/((1-rho^2)^2), 2*rho/((1-rho^2)^2)
    X = stats::toeplitz(c(4*rho,
                   -(rho^2 + 1),
                   rep(0,n - 2)))
    X[1,1] <-  2*rho
    X[n,n] <-  2*rho
    return((1/(1 - rho^2)^2) * X)
  }
  X =  stats::toeplitz( c(1 + rho^2, -rho, rep(0, n - 2)))
  X[1,1] <- X[n,n] <- 1
  return((1/(1 - rho^2)) * X)

}

#' Inverse of Compound Symmetric Matrix.
#'
#' @param n dimension of the matrix
#' @param rho compound symmetric factor
#' @param deriv logical, whether to return the derivative of the inverse.
#'
#' @return inverse matrix or derivative of the inverse matrix.
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' invCS(3,.25,TRUE)
#' }
invCS <- function(n, rho, deriv = FALSE){
  if (!(n > 1)) stop("n needs to be greater than 1")
  if (!(rho < 1 && rho > -1)) stop("rho needs to be < 1")
  alpha = sqrt(1 - rho)

  if (deriv) {
    X = diag(1/alpha^4, n) -
      matrix( ((n - 1) * rho^2 + 1)/((1 - rho)^2 * ((n - 1) * rho + 1)^2),
              nrow = n, ncol = n)
    return(X)
  }

  X = diag(1/alpha^2,n) - matrix(rho/(alpha^2 * (alpha^2 + n * rho)),
                                 nrow = n, ncol = n)
  return(X)
}



#' Positive Matrix Square Root
#' @description Computes a positive symmetric square root  matrix for a
#'    positive definite input matrix. Used in the inverted matrix variate
#'    t-distribution.
#' @param A positive definite p x p real-valued matrix.
#'
#' @return a symmetric square root matrix for A. ie, \eqn{B = t(B)} and
#'    \eqn{B \%*\% B = A}.
#' @export
#' @seealso \code{\link{posmatsqrtinv}}
#'
#' @examples
#' A = diag(5) + 1
#' B = posmatsqrt(A)
#' sum(abs(B - t(B)))
#' sum(abs(A - B %*% B))
posmatsqrt <- function(A) {
  # this isn't the fastest way and if you have to do this a lot find a
  # better way returns symmetric square root of A if it exists: B %*% B = A
  # does not test if A is positive definite
  if (!(is.numeric(A))) stop("Non-numeric input.")
  A <- as.matrix(A)
  if (!(dim(A)[1] == dim(A)[2]))
    stop("Matrix must be square. Dimensions: ", dim(A))

  e <- eigen(A, symmetric = TRUE)
  V <- e$vectors
  if (!(all(e$values > 0))) stop("Not all eigenvalues positive. e =",e$values)
  if (dim(A)[1] == 1) B <- (sqrt(e$values)) * V %*% t(V) else
    B <- V %*% diag(sqrt(e$values)) %*% t(V)
  return(B)
}

#' Positive Matrix Square Root Inverse
#' @description Computes the inverse of a positive symmetric square root  matrix for a
#'    positive definite input matrix. Used in the inverted matrix variate
#'    t-distribution.
#' @param A positive definite p x p real-valued matrix.
#'
#' @return a symmetric square root matrix for A. ie, \eqn{B = t(B)} and
#'    \eqn{B \%*\% B = A}.
#' @export
#' @seealso \code{\link{posmatsqrt}}
#'
#' @examples
#' A = diag(5) + 1
#' B = posmatsqrt(A)
#' C = posmatsqrtinv(A)
#' sum(abs(B - t(B)))
#' sum(abs(A - B %*% B))
#' B %*% C
posmatsqrtinv <- function(A) {
  # this isn't the fastest way and if you have to do this a lot find a
  # better way returns symmetric square root of A if it exists: B %*% B = A
  # does not test if A is positive definite
  if (!(is.numeric(A))) stop("Non-numeric input.")
  A <- as.matrix(A)
  if (!(dim(A)[1] == dim(A)[2]))
    stop("Matrix must be square. Dimensions: ", dim(A))

  e <- eigen(A, symmetric = TRUE)
  V <- e$vectors
  if (!(all(e$values > 0))) stop("Not all eigenvalues positive. e =",e$values)
  if (dim(A)[1] == 1) B <- (1/sqrt(e$values)) * V %*% t(V) else
    B <- V %*% diag(1/sqrt(e$values)) %*% t(V)
  return(B)
}

#' Cholesky of Random Wishart Distributed Matrices
#'
#' @description Generate n random matrices, distributed according to the Cholesky
#'     decomposition of a Wishart distribution with parameters \code{Sigma} and
#'     \code{df}, \eqn{W_p(Sigma, df)}.
#'
#' @param n integer sample size.
#' @param df numeric parameter, "degrees of freedom".
#' @param Sigma positive definite \eqn{(p * p)} "scale" matrix, the matrix parameter of the distribution.
#'
#' @return a numeric array, say R, of dimension \eqn{p * p * n}, where each \code{R[,,i]} is a Cholesky decomposition of a realization of the Wishart distribution \eqn{W_p(Sigma, df)}. Based on a modification of the existing code for the \code{rWishart} function
#'
#' @seealso \code{\link{rWishart}}, \code{\link{rInvCholWishart}}
#' @useDynLib matrixdist
#' @export
#'
#' @examples
#' # How it is parameterized:
#' set.seed(20180211)
#' A <- rCholWishart(1,10,3*diag(5))[,,1]
#'
#' set.seed(20180211)
#' B <- rInvCholWishart(1,10,3*diag(5))[,,1]
#'
#' A %*% B
#'
#' set.seed(20180211)
#' C <- chol(rWishart(1,10,3*diag(5))[,,1])
#'
rCholWishart <- function(n, df, Sigma){
  if (!is.numeric(Sigma))
    stop("Sigma must be numeric.")
  Sigma <- as.matrix(Sigma)
  dims = dim(Sigma)
  if (n < 1 || !(is.numeric(n)))
    stop("N must be 1 or larger.")
  if (!is.matrix(Sigma)  || dims[1] != dims[2])
    stop("'Sigma' must be a square, real matrix");
  if (!is.numeric(df) || df < dims[1] || dims[1] <= 0 )
    stop("inconsistent degrees of freedom and dimension")
  .Call("rCholWishart", n, df, Sigma, PACKAGE = "matrixdist")
}



#' Inverse of Cholesky of Random Wishart Distributed Matrices
#'
#' @description Generate n random matrices, distributed according to the Inverse of  Cholesky
#'     decomposition of a Wishart distribution with parameters \code{Sigma} and
#'     \code{df}, \eqn{W_p(Sigma, df)}.
#'
#' @param n integer sample size.
#' @param df numeric parameter, "degrees of freedom".
#' @param Sigma positive definite \eqn{(p * p)} "scale" matrix, the matrix parameter of the distribution.
#'
#' @return a numeric array, say R, of dimension \eqn{p * p * n}, where each \code{R[,,i]} is a Cholesky decomposition of a realization of the Wishart distribution \eqn{W_p(Sigma, df)}. Based on a modification of the existing code for the \code{rWishart} function
#'
#' @seealso \code{\link{rWishart}} and \code{\link{rCholWishart}}
#' @useDynLib matrixdist
#' @export
#'
#' @examples
#'
#' rInvCholWishart(1,10,diag(5))
#'
rInvCholWishart <- function(n, df, Sigma){
  if (!is.numeric(Sigma))
    stop("Sigma must be numeric.")
  Sigma <- as.matrix(Sigma)
  dims = dim(Sigma)
  if (n < 1 || !(is.numeric(n)))
    stop("N must be 1 or larger.")
  if (!is.matrix(Sigma)  || dims[1] != dims[2])
    stop("'Sigma' must be a square, real matrix");
  if (!is.numeric(df) || df < dims[1] || dims[1] <= 0 )
    stop("inconsistent degrees of freedom and dimension")
  .Call("rInvCholWishart", n, df, Sigma, PACKAGE = "matrixdist")
}



.onUnload <- function(libpath) {
  library.dynam.unload("matrixdist", libpath)
}
