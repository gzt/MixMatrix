#' matrixdist: a package for some random matrix distributions.
#'
#' The matrixdist package current provides for generating
#' matrix variate random variables, computing matrix variate densities,
#' and fitting the MLEs for parameters of matrix variate normal distributions.
#' Currently, normal, t, and inverted t distributions are supported for generation
#' and density computation.
#'
#' @docType package
#' @name matrixdist
NULL

#' rmatrixnorm: matrix variate normal
#'
#' @description  generate random draws from a matrix variate normal distribution
#'
#' @family matrixnorm
#' @param n number of observations to generate - must be a positive integer.
#' @param mean \eqn{p X q}  matrix of means
#' @param L \eqn{p X p}  matrix specifying relations among the rows.
#'    By default, an identity matrix.
#' @param R \eqn{q X q}  matrix specifying relations among the columns.
#'    By default, an identity matrix.
#' @param U \eqn{LL^T}  - \eqn{p X p}  positive definite variance-covariance
#'    matrix for rows, computed from \eqn{L} if not specified.
#' @param V \eqn{R^T R}  - \eqn{q X q}  positive definite variance-covariance
#'    matrix for columns, computed from \eqn{R}  if not specified.
#' @param list Defaults to \code{FALSE} . If this is \code{TRUE} , then the
#'    output will be a list of matrices.
#' @param array If \eqn{n = 1}  and this is not specified and \code{list}  is
#'    \code{FALSE} , the function will return a matrix containing the one
#'    observation. If \eqn{n > 1} , should be the opposite of \code{list} .
#'    If \code{list}  is \code{TRUE} , this will be ignored.
#' @return This returns either a list of \eqn{n}  \eqn{p X q}  matrices or
#'    a \eqn{p X q X n}  array.
#' @export
#'
#' @examples
#' set.seed(20180202)
#' rmatrixnorm(n=1,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'    L=matrix(c(2,1,0,.1),nrow=2),list=FALSE)
#' set.seed(20180202)
#' A <- rmatrixnorm(n=10,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'    L=matrix(c(2,1,0,.1),nrow=2),list=TRUE)
#' A[[1]]
#' set.seed(20180202)
#' B <- rmatrixnorm(n=10,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'    L=matrix(c(2,1,0,.1),nrow=2),list=FALSE)
#' B[ , , 1]
#'
rmatrixnorm <- function(n, mean, L = diag(dim(as.matrix(mean))[1]),
                        R = diag(dim(as.matrix(mean))[2]), U = L %*% t(L),
                        V = t(R) %*% R, list = FALSE, array = NULL) {
  if (!is.numeric(n)) stop("n is not numeric")
  if (!(n > 0)) stop("n must be > 0. n = ", n)
  mean <- as.matrix(mean)
  U <- as.matrix(U)
  V <- as.matrix(V)
  dims <- dim(mean)
  if (!(all(is.numeric(mean), is.numeric(U),is.numeric(V))))
    stop("Non-numeric input. ")

  # checks for conformable matrix dimensions
  if (!(dims[1] == dim(U)[2] && dim(U)[1] == dim(U)[2] &&
        dims[2] == dim(V)[1] && dim(V)[1] == dim(V)[2])) {
    stop("Non-conforming dimensions.", dims, dim(U),dim(V))
  }

  cholU <- chol(U)
  cholV <- chol(V)
  nobs = prod(dims)*n
  mat <- array(stats::rnorm(nobs), dim = c(dims,n))

  result <- array(apply(mat, 3, function(x) mean + t(cholU) %*% x %*% (cholV)),
                  dim = c(dims,n))
  if (n == 1 && list == FALSE && is.null(array)) {
    return(result[ , , 1])
    # if n = 1 and you don't specify arguments, if just returns a matrix
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


#' dmatrixnorm: matrix variate normal
#'
#' @description compute densities from a matrix variate normal distribution
#'
#' @family matrixnorm
#' @param x \eqn{p X q} input matrix
#' @param mean \eqn{p X q} matrix of means. By default, a matrix of \eqn{0}s
#'    with size taken from \code{x}
#' @param L \eqn{p X p} matrix specifying relations among the rows. By default,
#'    an identity matrix.
#' @param R \eqn{q X q} matrix specifying relations among the columns. By
#'    default, an identity matrix.
#' @param U \eqn{LL^T} - \eqn{p X p} positive definite variance-covariance matrix
#'    for rows, computed from \eqn{L} if not specified.
#' @param V \eqn{R^TR} - \eqn{q X q} positive definite variance-covariance matrix
#'    for columns, computed from \eqn{R} if not specified.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#'
#' @return Returns the density at the provided observation.
#' @export
#'
#' @examples
#' set.seed(20180202)
#' A <- rmatrixnorm(n=1,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#' L <- matrix(c(2,1,0,.1),nrow=2))
#' dmatrixnorm(A,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'   L=matrix(c(2,1,0,.1),nrow=2),log=TRUE )
#'
#'
dmatrixnorm <- function(x, mean = array(0L, dim(as.matrix(x))[1:2]),
                        L = diag(dim(mean)[1]),
                        R = diag(dim(mean)[2]), U = L %*% t(L),
                        V = t(R) %*% R, log = FALSE) {
  if (!(all(is.numeric(x), is.numeric(mean),
           is.numeric(U),is.numeric(V)))) stop("Non-numeric input. ")
    x <- as.matrix(x)
    mean <- as.matrix(mean)
    U <- as.matrix(U)
    V <- as.matrix(V)
    dims <- dim(x)
    if (!(dims[1] == dim(U)[2] && dim(U)[1] == dim(U)[2] &&
          dims[2] == dim(V)[1] && dim(V)[1] == dim(V)[2])) {
      stop("Non-conforming dimensions.", dims, dim(U),dim(V))
    }
    # you should be using small enough matrices that determinants and
    # inverses aren't a pain.  also presumes not using a singular matrix
    # normal distribution
    p <- dim(U)[1]  #square matrices so only need first dimension
    n <- dim(V)[1]
    detU <- det(U)
    detV <- det(V)
    if (!(detU > 0 && detV > 0)) stop("non-invertible matrix", detU, detV)
    Uinv <- solve(U)
    Vinv <- solve(V)
    XM <- x - mean
    logresult <- -0.5 * n * p * log(2 * pi) - 0.5 * n * log(detU) -
      0.5 * p * log(detV) - 0.5 * sum(diag(Vinv %*% t(XM) %*% Uinv %*% (XM)))
    if (log) {
        return(logresult)
      } else {
        return(exp(logresult))
      }
}

#' dmatrixnorm.test
#' @description Equivalent to dmatrixnorm except it works by unrolling
#'     to a vector.
#' @param x \eqn{p X q} input matrix
#' @param mean \eqn{p X q} matrix of means. By default, a matrix of \eqn{0}s
#'     with size taken from \code{x}
#' @param L \eqn{p X p} matrix specifying relations among the rows. By default,
#'     an identity matrix.
#' @param R \eqn{q X q} matrix specifying relations among the columns. By
#'    default, an identity matrix.
#' @param U \eqn{LL^T} - \eqn{p X p} positive definite variance-covariance
#'    matrix for rows, computed from \eqn{L} if not specified.
#' @param V \eqn{R^T R} - \eqn{q X q} positive definite variance-covariance
#'    matrix for columns, computed from \eqn{R} if not specified.
#' @param log Whether to return the density on the log scale.
#'
#' @return Returns the density at the provided observation. This is an
#'    alternative method of computing which works by flattening out into
#'    a vector instead of a matrix.
#' @keywords internal
#' #'
#' @examples
#' set.seed(20180202)
#' A <- rmatrixnorm(n=1,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'     L=matrix(c(2,1,0,.1),nrow=2))
#' \dontrun{dmatrixnorm.test (A,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'     L=matrix(c(2,1,0,.1),nrow=2),log=TRUE )}

dmatrixnorm.test <- function(x, mean = array(0L, dim(as.matrix(x))),
                             L = diag(dim(mean)[1]), R = diag(dim(mean)[2]),
                             U = L %*% t(L), V = t(R) %*% R, log = FALSE) {
    # results should equal other option - works by unrolling into MVN
  if (!(all(is.numeric(x), is.numeric(mean), is.numeric(L), is.numeric(R),
           is.numeric(U),is.numeric(V)))) stop("Non-numeric input. ")
    x <- as.matrix(x)
    mean <- as.matrix(mean)
    U <- as.matrix(U)
    V <- as.matrix(V)
    dims <- dim(x)
    if (!(dims[1] == dim(U)[2] && dim(U)[1] == dim(U)[2] &&
          dims[2] == dim(V)[1] && dim(V)[1] == dim(V)[2])) {
      stop("Non-conforming dimensions. ", dims, dim(U),dim(V))
    }
    vecx <- as.vector(x)
    meanx <- as.vector(mean)
    VU <- V %x% U
    # you should be using small enough matrices that determinants aren't a pain
    # also presumes not using a singular matrix normal distribution
    p <- dim(U)[1]  #square matrices so only need first dimension
    n <- dim(V)[1]
    detVU <- det(VU)
    if (!(detVU > 0)) stop("non-invertible matrix")
    UVinv <- solve(VU)
    XM <- vecx - meanx
    logresult <- -0.5 * n * p * log(2 * pi) - 0.5 * log(det(VU)) -
      0.5 * sum(diag(t(XM) %*% UVinv %*% XM))
    if (log) {
        return(logresult)
      } else {
        return(exp(logresult))
      }
}

#' mle.matrixnorm:
#'
#' @description Maximum likelihood estimation for matrix normal distributions
#'
#' @param data Either a list of matrices or a 3-D array with matrices in
#'    dimensions 2 and 3, indexed by dimension 1.
#' @param row.mean By default, \code{FALSE}. If \code{TRUE}, will fit a
#'    common mean within each row. If both this and \code{col.mean} are
#'    \code{TRUE}, there will be a common mean for the entire matrix.
#' @param col.mean By default, \code{FALSE}. If \code{TRUE}, will fit a
#'    common mean within each row. If both this and \code{row.mean} are
#'    \code{TRUE}, there will be a common mean for the entire matrix.
#' @param row.variance Imposes a variance structure on the rows. Either
#'     'none' or 'AR(1)'. Only positive correlations are allowed for AR(1).
#' @param col.variance  Imposes a variance structure on the columns.
#'     Either 'none' or 'AR(1)'. Only positive correlations are allowed for
#'     AR(1).
#' @param tol Convergence criterion. Measured against square deviation
#'    between iterations of the two variance-covariance matrices.
#' @param max.iter Maximum possible iterations of the algorithm.
#' @param U (optional) Can provide a starting point for the U matrix.
#'    By default, an identity matrix.
#' @param V (optional) Can provide a starting point for the V matrix.
#'    By default, an identity matrix.
#'
#' @return Returns a list with a mean matrix, a \eqn{U} matrix, a \eqn{V}
#'    matrix, the number of iterations, and error at the time of stopping.
#' @export
#'
#' @examples
#' set.seed(20180202)
#' A <- rmatrixnorm(n=100,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'    L=matrix(c(2,1,0,.1),nrow=2),list=TRUE)
#' results=mle.matrixnorm(A)
#' print(results)
#'
#'
mle.matrixnorm <- function(data, row.mean = FALSE, col.mean = FALSE,
                           row.variance = "none", col.variance = "none",
                           tol = 1e-09, max.iter = 100, U, V) {
    if (class(data) == "list") data <- array(unlist(data),
                                       dim = c(nrow(data[[1]]),
                                               ncol(data[[1]]), length(data)))
    if (!all(is.numeric(data),is.numeric(tol),
            is.numeric(max.iter))) stop("Non-numeric input. ")
    if (!(missing(U))){
      if (!(is.numeric(U))) stop("Non-numeric input.")
    }
    if (!(missing(V))){
      if (!(is.numeric(V))) stop("Non-numeric input.")
    }


    # if data is array, presumes indexed over first column (same as output
    # of rmatrixnorm) if list, presumes is a list of the matrices
    dims <- dim(data)

    if (max(dims[1]/dims[2], dims[2]/dims[1]) > (dims[3] - 1))
        warning("Need more observations to estimate parameters.")
    # don't have initial starting point for U and V, start with diag.
    if (missing(U))
        U <- diag(dims[1])
    if (missing(V))
        V <- diag(dims[2])
    mu <- apply(data, c(1, 2), mean)
    if (row.mean) {
        # should make it so that the mean is constant within a row
        mu <- matrix(apply(mu, 1, mean), nrow = dims[1], ncol = dims[2])
    }
    if (col.mean) {
        # should make it so that the mean is constant within a column
        mu <- matrix(apply(mu, 2, mean), nrow = dims[1], ncol = dims[2], byrow = T)
    }
    # if both are true, this should make it so the mean is constant all over
    swept.data <- sweep(data, c(1, 2), mu)
    iter <- 0
    error.term <- 1e+40
    if (col.variance == "AR(1)"){
      if (V[1,2] > 0) rho.col = V[1,2] else rho.col <- 0.3
    }
    if (row.variance == "AR(1)"){
      if (V[1,2] > 0) rho.row = V[1,2] else rho.row <- 0.3
    }

    while (iter < max.iter && error.term > tol) {
        # make intermediate matrix, then collapse to final version
        if (col.variance == "AR(1)") {
            var <- V[1, 1]
            nLL <- function(theta) {
              Vmat <- theta[1] * toepgenerate(dims[2], theta[1]) - sum(apply(swept.data, 3,
                  function(x) dmatrixnorm(x, log = TRUE, U = U, V = Vmat)))
            }
            fit0 <- stats::optim(c(var, rho.col), nLL, method = "L-BFGS-B",
                 hessian = FALSE, lower = c(0, 0), upper = c(Inf, 0.999))
            var <- fit0$par[1]
            rho.col <- fit0$par[2]
            new.V <- var * toepgenerate(dims[2], rho.col)
        } else {
            inter.V <- apply(swept.data, 3, function(x) (t(x) %*% solve(U) %*% x))
            # collapsed into a (row*column) * n, which is then gathered and fixed.
            new.V <- matrix(apply(inter.V, 1, sum),
                            nrow = dims[2])/(dims[3] * dims[1])
        }

        if (row.variance == "AR(1)") {
            nLL <- function(theta){
              Umat <- toepgenerate(dims[1],theta)
              -sum(apply(swept.data, 3,
                  function(x) dmatrixnorm(x, log = TRUE, V = new.V, U = Umat)))
            }
            fit0 <- stats::optim(rho.row, nLL, method = "L-BFGS-B",
                                hessian = FALSE, lower = 0, upper = 0.999)
            rho.row <- fit0$par
            new.U <- toepgenerate(dims[1], rho.row)
        } else {
            inter.U <- apply(swept.data, 3, function(x) ((x) %*% solve(V) %*% t(x)))
            # collapsed into a (row*column) * n, which is then gathered and fixed.
            new.U <- matrix(apply(inter.U, 1, sum), nrow = dims[1])/(dims[3] * dims[2])
            new.U <- new.U/(new.U[1, 1])
        }
        # only identifiable up to a constant, so have to fix something at 1
        # should perhaps change - makes doing other restrictions on variance
        # harder.  compute differences with prior iterations:
        error.term <- sum((new.V - V)^2) + sum((new.U - U)^2)
        V <- new.V
        U <- new.U
        iter <- iter + 1
    }
    if (iter >= max.iter || error.term > tol)
        warning("Failed to converge")

    return(list(mean = mu, U = U, V = V, iter = iter, tol = error.term))
}

#' toepgenerate
#'
#' @description generate AR(1) correlation matrices
#'
#' @param n number of columns/rows
#' @param rho correlation parameter
#'
#' @return Toeplitz $n x n$ matrix with 1 on the diagonal and $rho^k$ on
#'    the other diagonals, where $k$ is distance from the main diagonal.
#'    Used internally but it is useful for generating your own random matrices.
#' @keywords internal
#' @export
#'
#' @examples
#' rho <- .9
#' n <- 6
#' toepgenerate(n,rho)
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

