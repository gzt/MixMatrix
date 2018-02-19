

#' matrixnorm: matrix variate Normal distribution functions
#'
#' @description  Density and random generation for the matrix variate normal distribution
#'
#' @family matrixnorm
#' @param n number of observations to generate - must be a positive integer.
#' @param x quantile for density
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
#' @param force If TRUE, will take the input of \code{L} and/or \code{R}
#'    directly - otherwise computes \code{U} and \code{V} and uses Cholesky
#'    decompositions. Useful for generating degenerate normal distributions.
#'    Will also override concerns about potentially singular matrices
#'    unless they are not, in fact, invertible.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @return This returns either a list of \eqn{n}  \eqn{p X q}  matrices or
#'    a \eqn{p X q X n}  array.
#' @export
#'
#' @seealso \code{rnorm} \code{Distributions}
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
#' dmatrixnorm(A[[1]],mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'   L=matrix(c(2,1,0,.1),nrow=2),log=TRUE )
#'
rmatrixnorm <- function(n, mean,
                        L = diag(dim(as.matrix(mean))[1]),
                        R = diag(dim(as.matrix(mean))[2]),
                        U = L %*% t(L),
                        V = t(R) %*% R,
                        list = FALSE,
                        array = NULL,
                        force = FALSE) {
  if (!is.numeric(n)) stop("n is not numeric")
  if (!(n > 0)) stop("n must be > 0. n = ", n)
  mean <- as.matrix(mean)
  U <- as.matrix(U)
  V <- as.matrix(V)
  if (missing(L))
    if (!symm.check(U)) stop("U not symmetric.")
  if (missing(R))
    if (!symm.check(V)) stop("V not symmetric.")
  dims <- dim(mean)
  if (!(all(is.numeric(mean), is.numeric(U),is.numeric(V))))
    stop("Non-numeric input. ")

  # checks for conformable matrix dimensions
  if (!(dims[1] == dim(U)[2] && dim(U)[1] == dim(U)[2] &&
        dims[2] == dim(V)[1] && dim(V)[1] == dim(V)[2])) {
    stop("Non-conforming dimensions.", dims, dim(U),dim(V))
  }
  if (force && !missing(L)) cholU <- L else cholU <- chol.default(U)
  if (force && !missing(R)) cholV <- R else cholV <- chol.default(V)

  if (!force && (min(diag(cholU)) < 1e-6 || min(diag(cholV)) < 1e-6) ) {
      stop("Potentially singular covariance, use force = TRUE if intended. ",
           min(diag(cholU)), min(diag(cholV)))
  }
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


#' @describeIn rmatrixnorm Density calculation for matrix variate normal distributions.
#' @export
#'
dmatrixnorm <- function(x, mean = array(0, dim(as.matrix(x))[1:2]),
                        L = diag(dim(mean)[1]),
                        R = diag(dim(mean)[2]), U = L %*% t(L),
                        V = t(R) %*% R, log = FALSE) {
  if (!(all(is.numeric(x), is.numeric(mean),
           is.numeric(U),is.numeric(V)))) stop("Non-numeric input. ")
    x <- as.matrix(x)
    mean <- as.matrix(mean)
    U <- as.matrix(U)
    V <- as.matrix(V)
    if (!symm.check(U)) stop("U not symmetric.")
    if (!symm.check(V)) stop("V not symmetric.")
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
    cholU = chol.default(U)
    cholV = chol.default(V)
    detU <- prod(diag(cholU))^2
    detV <- prod(diag(cholV))^2
    if (!(detU > 1e-8 && detV > 1e-8)) stop("non-invertible matrix", detU, detV)
    Uinv <- chol2inv(cholU)
    Vinv <- chol2inv(cholV)
    XM <- x - mean
    logresult <- -0.5 * n * p * log(2 * pi) - 0.5 * n * log(detU) -
      0.5 * p * log(detV) - 0.5 * sum(diag( tcrossprod(Vinv, XM) %*% crossprod(Uinv, XM)))
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
#' @param unrolled logical, \code{FALSE} by default. If \code{x}
#'    is already unrolled, select \code{TRUE}. This will take the dimensions
#'    from the variance matrices, so they must be specified.
#' @param log logical - whether to return the density on the log scale.
#'
#' @return Returns the density at the provided observation. This is an
#'    alternative method of computing which works by flattening out into
#'    a vector instead of a matrix.
#'
#' @export
#'
#'
#' @examples
#' set.seed(20180202)
#' A <- rmatrixnorm(n=1,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'     L=matrix(c(2,1,0,.1),nrow=2))
#' \dontrun{dmatrixnorm.unroll (A,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'     L=matrix(c(2,1,0,.1),nrow=2),log=TRUE )}

dmatrixnorm.unroll <- function(x, mean = array(0L, dim(as.matrix(x))),
                             L = diag(dim(mean)[1]), R = diag(dim(mean)[2]),
                             U = L %*% t(L), V = t(R) %*% R, log = FALSE,
                             unrolled = FALSE) {
    # results should equal other option - works by unrolling into MVN
  if (!(all(is.numeric(x), is.numeric(mean), is.numeric(L), is.numeric(R),
           is.numeric(U),is.numeric(V)))) stop("Non-numeric input. ")
    x <- as.matrix(x)
    mean <- as.matrix(mean)
    U <- as.matrix(U)
    V <- as.matrix(V)
    if (!symm.check(U)) stop("U not symmetric.")
    if (!symm.check(V)) stop("V not symmetric.")
    dims <- dim(mean)
    if (unrolled) {
      dims <- c(dim(U)[1] ,dim(V)[1])
    }
    if (!(dims[1] == dim(U)[2] && dim(U)[1] == dim(U)[2] &&
          dims[2] == dim(V)[1] && dim(V)[1] == dim(V)[2])) {
      stop("Non-conforming dimensions. ", dims, dim(U),dim(V))
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
    p <- dim(U)[1]  #square matrices so only need first dimension
    n <- dim(V)[1]
    cholVU <- chol.default(VU)
    detVU <- prod(diag(cholVU))^2
    if (!(detVU > 1e-8)) stop("non-invertible matrix")
    UVinv <- chol2inv(cholVU)
    XM <- vecx - meanx
    logresult <- -0.5 * n * p * log(2 * pi) - 0.5 * log(detVU) -
      0.5 * sum(diag(t(XM) %*% UVinv %*% XM))
    if (log) {
        return(logresult)
      } else {
        return(exp(logresult))
      }
}

#' MLmatrixnorm:
#'
#' @description Maximum likelihood estimation for matrix normal distributions
#'
#' Maximum likelihood estimates exist for \eqn{N > max(p/q,q/p)+1} and are
#' unique for \eqn{N > max(p,q)}. This finds the estimate for the mean and then alternates
#' between estimates for the \eqn{U} and \eqn{V} matrices until convergence.
#' An AR(1) or compound symmetry restriction can be proposed for either or both
#' variance matrices. However, if they are inappropriate for the data, they may fail with
#' a warning.
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
#'     'none', 'AR(1)', or 'CS' for 'compound symmetry'.
#'     Only positive correlations are allowed for AR(1) and CS.
#'     Note that while maximum likelihood estimators are available (and used) for
#'     the unconstrained variance matrices, \code{optim} is used for any
#'     constraints so it will be considerably slower.
#' @param col.variance  Imposes a variance structure on the columns.
#'     Either 'none' or 'AR(1)' or 'CS'. Only positive correlations are allowed for
#'     AR(1) and CS.
#' @param tol Convergence criterion. Measured against square deviation
#'    between iterations of the two variance-covariance matrices.
#' @param max.iter Maximum possible iterations of the algorithm.
#' @param U (optional) Can provide a starting point for the U matrix.
#'    By default, an identity matrix.
#' @param V (optional) Can provide a starting point for the V matrix.
#'    By default, an identity matrix.
#' @param ... (optional) additional arguments can be passed to \code{optim}
#'    if using restrictions on the variance.
#'
#' @return Returns a list with a mean matrix, a \eqn{U} matrix, a \eqn{V}
#'    matrix, the variance parameter (the first entry of the variance matrices
#'    are constrained to be 1 for uniqueness), the number of iterations, the squared difference
#'    between iterations of the variance matrices at the time of stopping, the log likelihood,
#'    and a convergence code.
#' @export
#' @seealso \code{rmatrixnorm}
#'
#' @examples
#' set.seed(20180202)
#' A <- rmatrixnorm(n=100,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'    L=matrix(c(2,1,0,.1),nrow=2),list=TRUE)
#' results=MLmatrixnorm(A, tol = 1e-5)
#' print(results)
#'
#'
MLmatrixnorm <- function(data, row.mean = FALSE, col.mean = FALSE,
                           row.variance = "none", col.variance = "none",
                           tol = 10*.Machine$double.eps^0.5, max.iter = 100, U, V,...) {
    if (class(data) == "list") data <- array(unlist(data),
                                       dim = c(nrow(data[[1]]),
                                               ncol(data[[1]]), length(data)))
    if (!all(is.numeric(data),is.numeric(tol),
            is.numeric(max.iter))) stop("Non-numeric input. ")
    if (!(missing(U))) {
      if (!(is.numeric(U))) stop("Non-numeric input.")
    }
    if (!(missing(V))) {
      if (!(is.numeric(V))) stop("Non-numeric input.")
    }
    row.set.var = FALSE
    if (row.variance == "AR(1)" || row.variance == "CS" ) row.set.var = TRUE

    col.set.var = FALSE
    if (col.variance == "AR(1)" || col.variance == "CS" ) col.set.var = TRUE
    # if data is array, presumes indexed over third column (same as output
    # of rmatrixnorm) if list, presumes is a list of the matrices
    dims <- dim(data)

    if (max(dims[1]/dims[2], dims[2]/dims[1]) > (dims[3] - 1))
        warning("Need more observations to estimate parameters.")
    # don't have initial starting point for U and V, start with diag.
    if (missing(U))
        U <- diag(dims[1])
    if (missing(V))
        V <- diag(dims[2])
    # mu <- apply(data, c(1, 2), mean)
    mu <- rowMeans(data, dims = 2)
    if (row.mean) {
        # make it so that the mean is constant within a row
        mu <- matrix(rowMeans(mu), nrow = dims[1], ncol = dims[2])
    }
    if (col.mean) {
        # make it so that the mean is constant within a column
        mu <- matrix(colMeans(mu), nrow = dims[1], ncol = dims[2], byrow = T)
    }
    # if both are true, this makes it so the mean is constant all over
    swept.data <- sweep(data, c(1, 2), mu)
    iter <- 0
    error.term <- 1e+40
    if (col.set.var) {
      if (V[1,2] > 0) {
        rho.col <- V[1,2]
          } else {
            inter.V <- apply(swept.data, 3,
                        function(x) crossprod((x), chol2inv(chol.default(U))) %*% x)
            # collapsed into a (row*column) * n, which is then gathered and fixed.
            V <- matrix(rowSums(inter.V, dims = 1),
                            nrow = dims[2])/(dims[3] * dims[1])
            if (col.variance == "AR(1)") rho.col <- V[1,2]/V[1,1]
            if (col.variance == "CS") rho.col <- mean(V[1,]/V[1,1])
            if (rho.col > .9) rho.col <- .9
            if (rho.col < 0) rho.col <- 0
            V <- varmatgenerate(dims[2],rho.col,col.variance)
          }
    }

  if (row.set.var) {
      if (U[1,2] > 0) {
        rho.row <- U[1,2]
        } else {
          inter.U <- apply(swept.data, 3,
                        function(x) ((x) %*% tcrossprod(chol2inv(chol.default(V)),(x))))
          # collapsed into a (row*column) * n, which is then gathered and fixed.
          U <- matrix(rowSums(inter.U , dims = 1),
                      nrow = dims[1])/(dims[3] * dims[2])
          if (row.variance == "AR(1)") rho.row <- U[1,2]/U[1,1]
          if (row.variance == "CS") rho.row <- mean(U[1,]/U[1,1])
          if (rho.row > .9) rho.row <- .9
          if (rho.row < 0) rho.row = 0
          U <- varmatgenerate(dims[1],rho.row,row.variance)
        }
    }

    varflag = FALSE
    while (iter < max.iter && error.term > tol && (!varflag)) {

        # make intermediate matrix, then collapse to final version
        if (col.set.var) {
          var <- V[1,1]
          var <- sum(apply(matrix(swept.data,ncol = dims[3]),2,
                    function(x) crossprod(x,
                                chol2inv(chol.default(V/var)) %x% chol2inv(chol.default((U)))) %*% x)) / (prod(dims))
          # can write LL for AR(1) directly and differentiate
          # det(SIGMA) = (1-rho^2)^(n-1) -> d(det)/d\rho = -2*rho*(n-1)(1-rho^2)^(n-2)
          # inv(SIGMA) -> 1/(1-rho^2) * tridiag with -rho on off diags, 1+rho^2, 1 on main
          # derivs: -(rho^2+1)/(1-rho^2)^2  , 4*rho/((1-rho^2)^2), 2*rho/((1-rho^2)^2)
          tmp <- array(apply(swept.data, 3, function(x) t(x) %*% chol2inv(chol.default(U)) %*% (x) ),
                       dim = c(dims[2],dims[2],dims[3]))
          nLL <- function(theta) {
            Vmat <- varinv(dims[2],theta,TRUE, col.variance)/var # try it
            B <- matrix(rowSums(apply(tmp, 3,
                        function(x) Vmat %*% x )), nrow = dims[2])
            # solved derivative, need to find where this is zero:
            0.5 * dims[1] * dims[3] * vardet(dims[2], theta, TRUE, col.variance) -
              (.5 ) * sum(diag(B)) # problem was wrong constant

          }
          if (!isTRUE(sign(nLL(0)) * sign(nLL(.999)) <= 0)) {
            warning("Endpoints of derivative of likelihood do not have opposite sign. Check variance specification.")
            rho.col = 0
            varflag = TRUE
          } else {
          fit0 <- stats::uniroot(nLL, c(0,.999),...)
          rho.col <- fit0$root
          }
          new.V <- var * varmatgenerate(dims[2], rho.col,col.variance)
        } else {
            inter.V <- apply(swept.data, 3, function(x) (crossprod(x, chol2inv(chol.default(U))) %*% x))
            # collapsed into a (row*column) * n, which is then gathered and fixed.
            new.V <- matrix(apply(inter.V, 1, sum),
                            nrow = dims[2])/(dims[3] * dims[1])
        }

        if (row.set.var) {
          tmp <- array(apply(swept.data, 3, function(x) (x) %*% chol2inv(chol.default(new.V)) %*% t(x) ),
                       dim = c(dims[1],dims[1],dims[3]))
            nLL <- function(theta) {
              Umat <- varinv(dims[1],theta,TRUE, row.variance)
              B <- matrix(rowSums(apply(tmp, 3,
                          function(x) Umat %*% x)), nrow = dims[1])
              # solved derivative, need to find where this is zero:
              0.5 * dims[2] * dims[3] * vardet(dims[1], theta, TRUE, row.variance) -
                (.5 ) * sum(diag(B)) # problem was wrong constant
            }
            if (!isTRUE(sign(nLL(0)) * sign(nLL(.999)) <= 0)) {
              warning("Endpoints of derivative of likelihood do not have opposite sign. Check variance specification.")
              rho.row = 0
              varflag = TRUE
            } else {
            fit0 <- stats::uniroot(nLL, c(0,.999),...)
            rho.row <- fit0$root
            }
            new.U <- varmatgenerate(dims[1], rho.row,row.variance)
        } else {
            inter.U <- apply(swept.data, 3, function(x) (tcrossprod(x, chol2inv(chol.default(new.V))) %*%  t(x)))
            # collapsed into a (row*column) * n, which is then gathered and fixed.
            new.U <- matrix(rowSums(inter.U, dims = 1),
                            nrow = dims[1])/(dims[3] * dims[2])
            new.U <- new.U/(new.U[1, 1])
        }
        # only identifiable up to a constant, so have to fix something at 1
        # should perhaps change - makes doing other restrictions on variance
        # harder.  compute differences with prior iterations:
        error.term <- sum((new.V - V)^2) + sum((new.U - U)^2)
        V <- new.V
        U <- new.U

        # reset mu to account for mean estimation *if* restricted means
        # unsure about this - I think only an issue if *known* covar matrix
        # in which case estimation is biased unless you correct for it
        #
        if (col.mean || row.mean)
          mu <- rowMeans(data, dims=2)
        if (row.mean){
          invV <- chol2inv(chol.default(V))
          sumV <- sum(invV)
          mu <- matrix(mu %*% (rowSums(invV))/sumV, nrow = dims[1], ncol = dims[2])
          }
        if(col.mean){
          invU <- chol2inv(chol.default(U))
          sumU <- sum(invU)
          mu <- matrix(colSums(invU) %*% mu /sumU, nrow = dims[1], ncol = dims[2], byrow = T)
        }
        if (col.mean || row.mean)
          swept.data <- sweep(data, c(1, 2), mu)

        iter <- iter + 1
    }
    if (iter >= max.iter || error.term > tol || varflag)
        warning("Failed to converge")

    converged = !(iter >= max.iter || error.term > tol || varflag)
    logLik = 0
    for (i in seq(dims[3])) {
      logLik = logLik + dmatrixnorm(data[,,i], mu, U = U, V = V, log = TRUE)
    }
    return(list(mean = mu,
                U = U,
                V = V/V[1,1],
                var = V[1,1],
                iter = iter,
                tol = error.term,
                logLik = logLik,
                convergence = converged,
                call = match.call()))
}
