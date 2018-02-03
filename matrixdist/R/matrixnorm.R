#' matrixdist: a package for some random matrix distributions.
#'
#' The matrixdist package current provides for generating
#' matrix normal random variables, computing matrix normal densities,
#' and fitting the MLEs for parameters of matrix normal distributions.
#' In the future, this may expand to matrix t distributions or some
#' other matrix random variables, at least generation and distribution
#' if not parameter estimation.
#'
#' @docType package
#' @name matrixdist
NULL





#' rmatrixnorm.one
#' @param mean $p X q$ matrix of means
#' @param L $p X p$ matrix specifying relations among the rows. By default, an identity matrix.
#' @param R $q X q$ matrix specifying relations among the columns. By default, an identity matrix.
#' @param U $LL^T$ - $p X p$ positive definite variance-covariance matrix for rows, computed from $L$ if not specified.
#' @param V $R^TR$ - $q X q$ positive definite variance-covariance matrix for columns, computed from $R$ if not specified.
#'
#' @return Returns a matrix of one observation. This function is for internal use only.
#' @keywords internal
rmatrixnorm.one <- function(mean,L=diag(dim(mean)[1]),R=diag(dim(mean)[2]), U = L %*% t(L), V = t(R) %*% R){
  dims = dim(mean)
  # should probably do better error checking, checks for conformable matrix dimensions
  if(!(dims[1] == dim(U)[2] && dim(U)[1]==dim(U)[2] && dims[2]==dim(V)[1] && dim(V)[1]==dim(V)[2])) {
    stop("Non-conforming dimensions.")
  }
  n = prod(dims)
  mat = matrix(stats::rnorm(n),nrow=dims[1])
  cholU = chol(U)
  cholV = chol(V)
  result = mean + cholU %*% mat %*% t(cholV)
  dimnames(result) = dimnames(mean)
  return(result)
}



#' rmatrixnorm
#'
#' @param n number of observations to generate - must be a positive integer.
#' @param mean $p X q$ matrix of means
#' @param L $p X p$ matrix specifying relations among the rows. By default, an identity matrix.
#' @param R $q X q$ matrix specifying relations among the columns. By default, an identity matrix.
#' @param U $LL^T$ - $p X p$ positive definite variance-covariance matrix for rows, computed from $L$ if not specified.
#' @param V $R^TR$ - $q X q$ positive definite variance-covariance matrix for columns, computed from $R$ if not specified.
#' @param list Defaults to \code{FALSE}. If this is \code{TRUE}, then the output will be a list of matrices.
#' @param array If $n = 1$ and this is not specified and \code{list} is \code{FALSE}, the function will return a matrix containing the one observation. If $n > 1$, should be the opposite of \code{list}. If \code{list} is \code{TRUE}, this will be ignored.
#' @return This returns either a list of $n$ $p X q$ matrices or a $n X p X q$ array.
#' @export
#'
#' @examples
#'set.seed(20180202)
#'rmatrixnorm(n=1,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'    L=matrix(c(2,1,0,.1),nrow=2),list=FALSE)
#'set.seed(20180202)
#'A = rmatrixnorm(n=10,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'    L=matrix(c(2,1,0,.1),nrow=2),list=TRUE)
#'A[[1]]
#'set.seed(20180202)
#'B = rmatrixnorm(n=10,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'    L=matrix(c(2,1,0,.1),nrow=2),list=FALSE)
#'B[1 , , ]
#'
rmatrixnorm <- function(n, mean,L=diag(dim(mean)[1]),R=diag(dim(mean)[2]),
                        U = L %*% t(L), V = t(R) %*% R,list=FALSE,array=NULL ){
    if(!(n>0)) stop("n must be > 0")
  dims = dim(mean)
  if(n ==1 && list == FALSE && is.null(array)){
    return(rmatrixnorm.one(mean,U=U,V=V))
    #if n = 1 and you don't specify arguments, if just returns a matrix
  }
    if(list){
    return(lapply(1:n, function(x) rmatrixnorm.one(mean,U=U,V=V)))
    }
  if (!(list) && !(is.null(array))) {
    if(!(array)) warning("list FALSE and array FALSE, returning array")
  }
  result = array(data=NA,dim=c(n,dims), dimnames=list(NULL,dimnames(mean)))
  for(i in 1:n){
    # note this indexes by the first coord - use aperm() on results if you don't want that
    result[i,,] =  rmatrixnorm.one(mean,U=U,V=V)
  }
  return(result)
}


#' dmatrixnorm
#'
#' @param x $p X q$ input matrix
#' @param mean $p X q$ matrix of means. By default, a matrix of $0$s with size taken from \code{x}
#' @param L $p X p$ matrix specifying relations among the rows. By default, an identity matrix.
#' @param R $q X q$ matrix specifying relations among the columns. By default, an identity matrix.
#' @param U $LL^T$ - $p X p$ positive definite variance-covariance matrix for rows, computed from $L$ if not specified.
#' @param V $R^TR$ - $q X q$ positive definite variance-covariance matrix for columns, computed from $R$ if not specified.
#' @param log Whether to return the density on the log scale.
#'
#' @return Returns the density at the provided observation.
#' @export
#'
#' @examples
#' set.seed(20180202)
#' A = rmatrixnorm(n=1,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#' L=matrix(c(2,1,0,.1),nrow=2))
#' dmatrixnorm(A,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'   L=matrix(c(2,1,0,.1),nrow=2),log=TRUE )
#'
#'
dmatrixnorm <- function(x, mean = array(0L, dim(x)), L = diag(dim(x)[1]), R=diag(dim(x)[2]),
                        U = L %*% t(L), V = t(R) %*% R,log = FALSE){
  dims = dim(x)
  if(!(dims[1] == dim(U)[2] && dim(U)[1]==dim(U)[2] && dims[2]==dim(V)[1] && dim(V)[1]==dim(V)[2])) {
    stop("Non-conforming dimensions.")
  }
  # you should be using small enough matrices that determinants and inverses aren't a pain.
  # also presumes not using a singular matrix normal distribution
  p = dim(U)[1] #square matrices so only need first dimension
  n = dim(V)[1]
  detU = det(U)
  detV = det(V)
  if(!(detU>0||detV>0)) stop("non-invertible matrix")
  Uinv = solve(U)
  Vinv = solve(V)
  XM = x-mean
  logresult = -.5*n*p*log(2*pi) - .5*n*log(detU) - .5*p*log(detV) - .5*sum(diag(Vinv %*% t(XM) %*% Uinv %*% (XM)   ))
  if(log) return(logresult)
  else return(exp(logresult))
}

#' dmatrixnorm.test
#'
#' @param x $p X q$ input matrix
#' @param mean $p X q$ matrix of means. By default, a matrix of $0$s with size taken from \code{x}
#' @param L $p X p$ matrix specifying relations among the rows. By default, an identity matrix.
#' @param R $q X q$ matrix specifying relations among the columns. By default, an identity matrix.
#' @param U $LL^T$ - $p X p$ positive definite variance-covariance matrix for rows, computed from $L$ if not specified.
#' @param V $R^TR$ - $q X q$ positive definite variance-covariance matrix for columns, computed from $R$ if not specified.
#' @param log Whether to return the density on the log scale.
#'
#' @return Returns the density at the provided observation. This is an alternative method of computing which works by flattening out into a vector instead of a matrix.
#' @keywords internal
#' @examples
#' set.seed(20180202)
#' A = rmatrixnorm(n=1,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'     L=matrix(c(2,1,0,.1),nrow=2))
#' dmatrixnorm(A,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'     L=matrix(c(2,1,0,.1),nrow=2),log=TRUE )

dmatrixnorm.test <- function(x, mean = array(0L, dim(x)), L = diag(dim(x)[1]),
                    R=diag(dim(x)[2]),U = L %*% t(L), V = t(R) %*% R,log = FALSE){
  #results should equal other option - works by unrolling into MVN
  dims = dim(x)
  if(!(dims[1] == dim(U)[2] && dim(U)[1]==dim(U)[2] && dims[2]==dim(V)[1] && dim(V)[1]==dim(V)[2])) {
    stop("Non-conforming dimensions.")
  }
  vecx = as.vector(x)
  meanx = as.vector(mean)
  VU = V %x% U
  # you should be using small enough matrices that determinants aren't a pain.
  # also presumes not using a singular matrix normal distribution
  p = dim(U)[1] #square matrices so only need first dimension
  n = dim(V)[1]
  detVU = det(VU)
  if (!(detVU>0)) stop("non-invertible matrix")
  UVinv = solve(VU)
  XM = vecx-meanx
  logresult = -.5*n*p*log(2*pi) - .5*log(det(VU)) - .5*sum(diag( t(XM) %*% UVinv %*% XM   ))
  if(log) return(logresult)
  else return(exp(logresult))
}

#' mle.matrixnorm
#'
#' @param data Either a list of matrices or a 3-D array with matrices in dimensions 2 and 3, indexed by dimension 1.
#' @param row.restrict At present, no options available for this. In the future, intend to have options for, eg, AR(1) structure.
#' @param col.restrict  At present, no options available for this. In the future, intend to have options for, eg, AR(1) structure.
#' @param tol Convergence criterion. Measured against square deviation between iterations of the two variance-covariance matrices.
#' @param max.iter Maximum possible iterations of the algorithm.
#' @param U (optional) Can provide a starting point for the U matrix. By default, an identity matrix.
#' @param V (optional) Can provide a starting point for the V matrix. By default, an identity matrix.
#'
#' @return Returns a list with a mean matrix, a $U$ matrix, a $V$ matrix, the number of iterations, and error at the time of stopping.
#' @export
#'
#' @examples
#'set.seed(20180202)
#'A = rmatrixnorm(n=100,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'    L=matrix(c(2,1,0,.1),nrow=2),list=TRUE)
#' results=mle.matrixnorm(A)
#' print(results)
#'
#'
mle.matrixnorm = function(data,row.restrict="none",col.restrict="none",tol = 1e-9,max.iter=100,U,V){
  if(class(data) == "list") data = aperm(array(unlist(data), dim = c(nrow(data[[1]]), ncol(data[[1]]), length(data))),perm=c(3,1,2))

  # intend to implement toeplitz restriction later
  # if data is array, presumes indexed over first column (same as output of rmatrixnorm)
  # if list, presumes is a list of the matrices
  # will start by working with assumption that data is an array
  dims = dim(data)

  if(max(dims[2]/dims[3],dims[3]/dims[2])<(dims[1]-1)) warning("Need more observations to estimate parameters.")
  # don't have initial starting point for U and V, start with diag.
  if(missing(U)) U = diag(dims[2])
  if(missing(V)) V = diag(dims[3])
  mu = apply(data,c(2,3),mean)
  swept.data = sweep(data,c(2,3),mu)
  iter = 0
  error.term = 1e40

  while(iter < max.iter && error.term > tol){
    #make intermediate matrix, then collapse to final version
    inter.V = apply(swept.data, 1, function(x) (t(x) %*% solve(U) %*% x))
    new.V = matrix(apply( inter.V, 1, sum),nrow=dims[3]) / (dims[1]*dims[2])

    inter.U =  apply(swept.data, 1, function(x) ((x) %*% solve(V) %*% t(x)) )
    new.U = matrix(apply(inter.U ,1,sum),nrow=dims[2]) / (dims[1]*dims[3])
    new.U = new.U/(new.U[1,1]) # only identifiable up to a constant, so have to fix something at 1
    # compute differences with prior iterations:
    error.term = sum((new.V - V)^2) + sum((new.U - U)^2)
    V = new.V
    U = new.U
    iter = iter + 1
  }
  if(iter >= max.iter || error.term > tol) warning("Failed to converge")

  return(list(mean = mu, U=U, V = V, iter=iter, tol=error.term))
}

