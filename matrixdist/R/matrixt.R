

#' rmatrixt.one
#'
#' @param mean \eqn{p X q} This is really a "shift" rather than a mean as this
#'    is a central T, though the expected value will be equal to this if \eqn{df > 2}
#' @param df Degrees of freedom parameter
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
rmatrixt.one <- function(mean = matrix(0,nrow=2,ncol=2),df,L=diag(dim(mean)[1]),R=diag(dim(mean)[2]), U = L %*% t(L), V = t(R) %*% R){
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
  # is multiple U by the degrees of freedom.
  # somebody please correct me on this.
  USigma = stats::rWishart(1,df+dims[1]-1,solve(df*U))
  cholU = t(chol(solve(USigma[,,1])))
  cholV = chol(V)
  result = mean + t(cholU) %*% mat %*% (cholV)
  dimnames(result) = dimnames(mean)
  return(result)
}

#' Title
#'
#' @param n number of observations to generate
#' @param mean \eqn{p X q} This is really a "shift" rather than a mean as this
#'    is a central T, though the expected value will be equal to this if \eqn{df > 2}
#' @param df Degrees of freedom parameter
#' @param L \eqn{p X p}  matrix specifying relations among the rows. By default, an identity matrix.
#' @param R \eqn{q X q}  matrix specifying relations among the columns. By default, an identity matrix.
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
#' rmatrixt(n=1,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),df=10,
#'    L=matrix(c(2,1,0,.1),nrow=2),list=FALSE)
#' set.seed(20180202)
#' A = rmatrixt(n=10,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),df=10,
#'    L=matrix(c(2,1,0,.1),nrow=2),list=TRUE)
#' A[[1]]
#' set.seed(20180202)
#' B = rmatrixt(n=10,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),df=10,
#'    L=matrix(c(2,1,0,.1),nrow=2),list=FALSE)
#' B[1 , , ]
#' summary(rt(n=100,df=10))
#' summary(rmatrixt(n=100,matrix(0),df=10))
#'

rmatrixt <- function(n, mean,df,L=diag(dim(mean)[1]),R=diag(dim(mean)[2]),
                              U = L %*% t(L), V = t(R) %*% R,list=FALSE,array=NULL ){
  if(!(n>0)) stop("n must be > 0")
  dims = dim(mean)
  if(n ==1 && list == FALSE && is.null(array)){
    return(rmatrixt.one(mean,df,U=U,V=V))
    #if n = 1 and you don't specify arguments, if just returns a matrix
  }
  if(list){
    return(lapply(1:n, function(x) rmatrixt.one(mean,df,U=U,V=V)))
  }
  if (!(list) && !(is.null(array))) {
    if(!(array)) warning("list FALSE and array FALSE, returning array")
  }
  result = array(data=NA,dim=c(n,dims), dimnames=list(NULL,dimnames(mean)))
  for(i in 1:n){
    # note this indexes by the first coord - use aperm() on results if you don't want that
    result[i,,] =  rmatrixt.one(mean,df,U=U,V=V)
  }
  return(result)
}
