#' rmatrixt.one
#'
#' @param df  degrees of freedom (\eqn{>0}, maybe non-integer), \code{df = Inf} is allowed.
#' @param mean \eqn{p X q} This is really a "shift" rather than a mean as this
#'    is a central T, though the expected value will be equal to this if \eqn{df > 2}
#' @param L \eqn{p X p}  matrix specifying relations among the rows. By default, an identity matrix.
#' @param R \eqn{q X q}  matrix specifying relations among the columns. By default, an identity matrix.
#' @param U \eqn{LL^T}  - \eqn{p X p}  positive definite matrix for rows, computed
#'    from \eqn{L} if not specified. Note this is not the variance.
#' @param V \eqn{R^T R}  - \eqn{q X q}  positive definite matrix for columns,
#'    computed from \eqn{R}  if not specified.  Note this is not the variance.
#'
#' @return Returns a matrix of one observation. This function is for internal use only.
#' @keywords internal
#'
rmatrixt.one <- function(df, mean = matrix(0,nrow=2,ncol=2), L=diag(dim(mean)[1]),
                         R=diag(dim(mean)[2]), U = L %*% t(L), V = t(R) %*% R){
  mean = as.matrix(mean)
  U = as.matrix(U)
  V = as.matrix(V)
  dims = dim(mean)
  # should probably do better error checking, checks for conformable matrix dimensions
  if(!(dims[1] == dim(U)[2] && dim(U)[1]==dim(U)[2] && dims[2]==dim(V)[1] && dim(V)[1]==dim(V)[2])) {
    stop("Non-conforming dimensions.",dims,dim(U),dim(V))
  }
  n = prod(dims)
  mat = matrix(stats::rnorm(n),nrow=dims[1])
  # okay, here's the deal: the way this was formulated in the book
  # would have solve(U) here, but after working out the math it is clear
  # that what we have to do to get what is "usually" expected on this count
  # is multiply U by the degrees of freedom.
  # somebody please correct me on this if wrong.
  USigma = stats::rWishart(1,df+dims[1]-1,solve(df*U))[,,1]
  cholU = (chol(solve(USigma)))
  cholV = chol(V)
  result = mean + t(cholU) %*% mat %*% (cholV)
  dimnames(result) = dimnames(mean)
  return(result)
}

#' rmatrixt
#' @family matrixt
#' @description Random generation for the matrix variate t-distribution
#' @param n number of observations to generate
#' @param df  degrees of freedom (\eqn{>0}, maybe non-integer), \code{df = Inf} is allowed.
#' @param mean \eqn{p X q} This is really a "shift" rather than a mean as this
#'    is a central T, though the expected value will be equal to this if \eqn{df > 2}
#' @param L \eqn{p X p}  matrix specifying relations among the rows.
#'     By default, an identity matrix.
#' @param R \eqn{q X q}  matrix specifying relations among the columns.
#'     By default, an identity matrix.
#' @param U \eqn{LL^T}  - \eqn{p X p}  positive definite matrix for rows, computed
#'    from \eqn{L} if not specified. Note this is not the variance.
#' @param V \eqn{R^T R}  - \eqn{q X q}  positive definite matrix for columns,
#'    computed from \eqn{R}  if not specified.  Note this is not the variance.
#' @param list Defaults to \code{FALSE} . If this is \code{TRUE} , then the output will be a list of matrices.
#' @param array If \eqn{n = 1}  and this is not specified and \code{list}  is \code{FALSE} , the function will return a matrix containing the one observation. If \eqn{n > 1} , should be the opposite of \code{list} . If \code{list}  is \code{TRUE} , this will be ignored.
#' @return This returns either a list of \eqn{n}  \eqn{p X q}  matrices or a \eqn{n X p X q}  array.
#' @export
#'
#' @examples
#' set.seed(20180202)
#' rmatrixt(n=1,df=10,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'    L=matrix(c(2,1,0,.1),nrow=2),list=FALSE)
#' set.seed(20180202)
#' A = rmatrixt(n=10,df=10,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'    L=matrix(c(2,1,0,.1),nrow=2),list=TRUE)
#' A[[1]]
#' set.seed(20180202)
#' B = rmatrixt(n=10,df=10,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'    L=matrix(c(2,1,0,.1),nrow=2),list=FALSE)
#' B[1 , , ]
#' summary(rt(n=100,df=10))
#' summary(rmatrixt(n=100,df=10,matrix(0)))
#'

rmatrixt <- function(n, df, mean,L=diag(dim(mean)[1]),R=diag(dim(mean)[2]),
                              U = L %*% t(L), V = t(R) %*% R,list=FALSE,array=NULL ){
  if(!(n>0)) stop("n must be > 0")
  mean = as.matrix(mean)
  U = as.matrix(U)
  V = as.matrix(V)
  dims = dim(mean)
  if(n ==1 && list == FALSE && is.null(array)){
    return(rmatrixt.one(mean=mean,df=df,U=U,V=V))
    #if n = 1 and you don't specify arguments, if just returns a matrix
  }
  if(list){
    return(lapply(1:n, function(x) rmatrixt.one(mean=mean,df=df,U=U,V=V)))
  }
  if (!(list) && !(is.null(array))) {
    if(!(array)) warning("list FALSE and array FALSE, returning array")
  }
  result = array(data=NA,dim=c(n,dims), dimnames=list(NULL,dimnames(mean)))
  for(i in 1:n){
    # note this indexes by the first coord - use aperm() on results if you don't want that
    result[i,,] =  rmatrixt.one(mean=mean,df=df,U=U,V=V)
  }
  return(result)
}

#' dmatrixt
#'
#' @family matrixt
#' @description Density function for the matrix variate t-distribution
#' @param x \eqn{p X q} input matrix
#' @param df degrees of freedom (\eqn{>0}, maybe non-integer), \code{df = Inf} is allowed.
#' @param mean by default 0, matrix of same size as observations \code{x}
#' @param L \eqn{p X p}  matrix specifying relations among the rows. By default, an identity matrix.
#' @param R \eqn{q X q}  matrix specifying relations among the columns. By default, an identity matrix.
#' @param U \eqn{LL^T}  - \eqn{p X p}  positive definite matrix for rows, computed
#'    from \eqn{L} if not specified. Note this is not the variance.
#' @param V \eqn{R^T R}  - \eqn{q X q}  positive definite matrix for columns,
#'    computed from \eqn{R}  if not specified.  Note this is not the variance.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#'
#' @return returns the density.
#' @export
#'
#' @examples
#' set.seed(20180204)
#' x = rmatrixt(n=1,mean=matrix(0),df=1)
#' dt(x,1)
#' dmatrixt(x,df=1)
#'
dmatrixt <- function(x, df, mean = array(0L, dim(x)[1:2]), L = diag(dim(x)[1]),
                      R=diag(dim(x)[2]), U = L %*% t(L), V = t(R) %*% R,log = FALSE){
  x = as.matrix(x)
  mean = as.matrix(mean)
  U = as.matrix(U)
  V = as.matrix(V)
  dims = dim(x)
  if(!(dims[1] == dim(U)[2] && dim(U)[1]==dim(U)[2] && dims[2]==dim(V)[1] && dim(V)[1]==dim(V)[2])) {
    stop("Non-conforming dimensions.")
  }
  xm = x - mean
  # breaking equation into two parts:
  # the integrating constants (gammas) and the matrix algebra parts (mats)
  # done on the log scale
  #NB: I provided the correction to this that I did for rmatrixt as well (ie scale by df)
  gammas = lmvgamma((0.5)*(df+sum(dims)-1),dims[1]) -
            .5*prod(dims)*log(pi) - lmvgamma(.5*(df + dims[1]-1),dims[1])
  mats = -.5*dims[2]*log(det(df*U))-.5*dims[1]*log(det(V))-
          0.5*(df+sum(dims)-1)*log(det(diag(dims[1])+
          solve(df*U) %*% xm %*% solve(V) %*% t(xm)))
  results = gammas + mats
  if(log) return(results)
    else return(exp(results))
}

#' lmvgamma
#'
#' @description A special mathematical function related to the gamma function,
#'     generalized for multivariate gammas.
#'
#' @param x non-negative numeric vector, matrix, or array
#' @param p positive integer, dimension of a square matrix
#'
#' @return log of multivariate gamma for each entry of x. For non-log variant, see mvgamma.
#'
#' @seealso mvgamma
#' @export
#'
#' @examples
#' lgamma(1:12)
#' lmvgamma(1:12,1)
lmvgamma <- function(x,p){
  # p only makes sense as an integer but not testing that
  # x *could* be less than zero - same domain as gamma function
  # making sure that object returned is same shape as object passed
  dims = if(is.vector(x)) length(x) else dim(as.array(x))
  if(p<1) stop("p must be greater than or equal to 1. p = ",p)
  if(any(x <= 0)) stop("x must be greater than 0. x = ",x)
  i = seq_along(p)
  # the sum is why we have to do this awkwardly rather than just passing through
  result <- sapply(x, function(y) (p*(p-1)/4)*log(pi) + sum(lgamma(y + (1-i)/2)))
  return(array(result,dim=dims))
}

mvgamma <- function(x,p) exp(lmvgamma(x,p))

#' posmatsqrt
#' @description Computes a positive symmetric square root  matrix for a positive definite input matrix. Used in the inverted matrix variate t-distribution.
#' @param A positive definite p x p real-valued matrix.
#'
#' @return a symmetric square root matrix for A. ie, \eqn{B = t(B)} and \eqn{B \%*\% B = A}.
#' @export
#'
#' @examples
#' A = diag(5) + 1
#' B = posmatsqrt(A)
#' sum(abs(B - t(B)))
#' sum(abs(A - B %*% B))
posmatsqrt <- function(A){
  # this isn't the fastest way and if you have to do this a lot find a better way
  # returns symmetric square root of A if it exists: B %*% B = A
  # does not test if A is positive definite
  A = as.matrix(A)
  if(!(dim(A)[1]==dim(A)[2])) stop("Matrix must be square. Dimensions: ",dim(A))

  e <- eigen(A)
  V <- e$vectors
  B <- V %*% diag(sqrt(e$values)) %*% t(V)
  return(B)
}


#' rmatrixinvt.one
#'
#' @description Generate random draws from the inverted matrix variate t-distribution.
#' @param df  degrees of freedom (\eqn{>0}, maybe non-integer), \code{df = Inf} is allowed.
#' @param mean \eqn{p X q} This is really a "shift" rather than a mean as this
#'    is a central T, though the expected value will be equal to this if \eqn{df > 2}
#' @param L \eqn{p X p}  matrix specifying relations among the rows. By default, an identity matrix.
#' @param R \eqn{q X q}  matrix specifying relations among the columns. By default, an identity matrix.
#' @param U \eqn{LL^T}  - \eqn{p X p}  positive definite matrix for rows, computed
#'    from \eqn{L} if not specified. Note this is not the variance.
#' @param V \eqn{R^T R}  - \eqn{q X q}  positive definite matrix for columns,
#'    computed from \eqn{R}  if not specified.  Note this is not the variance.
#'
#' @return Returns a matrix of one observation. This function is for internal use only.
#' @keywords internal
#'
rmatrixinvt.one <- function(df, mean = matrix(0,nrow=2,ncol=2), L=diag(dim(mean)[1]),
                        R=diag(dim(mean)[2]), U = L %*% t(L), V = t(R) %*% R){
  mean = as.matrix(mean)
  U = as.matrix(U)
  V = as.matrix(V)
  dims = dim(mean)
  # should probably do better error checking, checks for conformable matrix dimensions
  if(!(dims[1] == dim(U)[2] && dim(U)[1]==dim(U)[2] && dims[2]==dim(V)[1] && dim(V)[1]==dim(V)[2])) {
    stop("Non-conforming dimensions.",dims,dim(U),dim(V))
  }
  n = prod(dims)
  mat = matrix(stats::rnorm(n),nrow=dims[1])
  S = stats::rWishart(1,df+dims[1]-1,diag(dims[1]))[,,1]
  Usqrt = posmatsqrt(U)
  Vsqrt = posmatsqrt(V)
  SXX = solve(posmatsqrt(S + mat %*% t(mat)))
  result = Usqrt %*% SXX %*% mat %*% Vsqrt + mean
  return(result)
}

#' rmatrixinvt
#'
#' @family matrixt
#' @description Generate random draws from the inverted matrix variate t-distribution
#' @param n number of observations
#' @param df  degrees of freedom (\eqn{>0}, maybe non-integer), \code{df = Inf} is allowed.
#' @param mean \eqn{p X q} This is really a "shift" rather than a mean.
#' @param L \eqn{p X p}  matrix specifying relations among the rows. By default, an identity matrix.
#' @param R \eqn{q X q}  matrix specifying relations among the columns. By default, an identity matrix.
#' @param U \eqn{LL^T}  - \eqn{p X p}  positive definite matrix for rows, computed
#'    from \eqn{L} if not specified. Note this is not the variance.
#' @param V \eqn{R^T R}  - \eqn{q X q}  positive definite matrix for columns,
#'    computed from \eqn{R}  if not specified.  Note this is not the variance.
#'
#' @return This returns either a list of \eqn{n}  \eqn{p X q}  matrices or a \eqn{n X p X q}  array.
#' @export
#'
rmatrixinvt <- function(n,df, mean = matrix(0,nrow=2,ncol=2), L=diag(dim(mean)[1]),
                            R=diag(dim(mean)[2]), U = L %*% t(L), V = t(R) %*% R){
  if(!(n>0)) stop("n must be > 0")
  mean = as.matrix(mean)
  U = as.matrix(U)
  V = as.matrix(V)
  dims = dim(mean)
  if(n ==1 && list == FALSE && is.null(array)){
    return(rmatrixinvt.one(mean,df,U=U,V=V))
    #if n = 1 and you don't specify arguments, if just returns a matrix
  }
  if(list){
    return(lapply(1:n, function(x) rmatrixinvt.one(mean,df,U=U,V=V)))
  }
  if (!(list) && !(is.null(array))) {
    if(!(array)) warning("list FALSE and array FALSE, returning array")
  }
  result = array(data=NA,dim=c(n,dims), dimnames=list(NULL,dimnames(mean)))
  for(i in 1:n){
    # note this indexes by the first coord - use aperm() on results if you don't want that
    result[i,,] =  rmatrixinvt.one(mean,df,U=U,V=V)
  }
  return(result)

}

#' dmatrixinvt
#' @family matrixt
#' @description Compute densities for the inverted matrix variate t-distribution.
#' @param x \eqn{p X q} input matrix
#' @param df degrees of freedom (\eqn{>0}, maybe non-integer), \code{df = Inf} is allowed.
#' @param mean by default 0, matrix of same size as observations \code{x}
#' @param L \eqn{p X p}  matrix specifying relations among the rows. By default, an identity matrix.
#' @param R \eqn{q X q}  matrix specifying relations among the columns. By default, an identity matrix.
#' @param U \eqn{LL^T}  - \eqn{p X p}  positive definite matrix for rows, computed
#'    from \eqn{L} if not specified. Note this is not the variance.
#' @param V \eqn{R^T R}  - \eqn{q X q}  positive definite matrix for columns,
#'    computed from \eqn{R}  if not specified.  Note this is not the variance.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#'
#' @return returns the density.
dmatrixinvt <- function(x, df, mean = array(0L, dim(x)[1:2]), L = diag(dim(x)[1]),
                        R=diag(dim(x)[2]), U = L %*% t(L), V = t(R) %*% R,log = FALSE){
  x = as.matrix(x)
  mean = as.matrix(mean)
  U = as.matrix(U)
  V = as.matrix(V)
  dims = dim(x)
  if(!(dims[1] == dim(U)[2] && dim(U)[1]==dim(U)[2] && dims[2]==dim(V)[1] && dim(V)[1]==dim(V)[2])) {
    stop("Non-conforming dimensions.")
  }
  xm = x - mean
  # breaking equation into two parts:
  # the integrating constants (gammas) and the matrix algebra parts (mats)
  # done on the log scale
  gammas = lmvgamma((0.5)*(df+sum(dims)-1),dims[1]) -
    .5*prod(dims)*log(pi) - lmvgamma(.5*(df + dims[1]-1),dims[1])
  mats = -.5*dims[2]*log(det(U))-.5*dims[1]*log(det(V))-
    0.5*(df-2)*log(det(diag(dims[1])-
                                   solve(U) %*% xm %*% solve(V) %*% t(xm)))
  results = gammas + mats
  if(log) return(results)
  else return(exp(results))
}
