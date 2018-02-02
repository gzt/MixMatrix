

rmatrixnorm.one <- function(mean,L=diag(dim(mean)[1]),R=diag(dim(mean)[2]), U = L %*% t(L), V = t(R) %*% R){
  dims = dim(mean)
  # should probably do better error checking, checks for conformable matrix dimensions
  if(!(dims[1] == dim(U)[2] && dim(U)[1]==dim(U)[2] && dims[2]==dim(V)[1] && dim(V)[1]==dim(V)[2])) {
    stop("Non-conforming dimensions.")
  }
  n = prod(dims)
  mat = matrix(rnorm(n),nrow=dims[1])
  cholU = chol(U)
  cholV = chol(V)
  result = mean + cholU %*% mat %*% t(cholV)
  dimnames(result) = dimnames(mean)
  return(result)
}



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


dmatrixnorm.test <- function(x, mean = array(0L, dim(x)), L = diag(dim(x)[1]), R=diag(dim(x)[2]), 
                        U = L %*% t(L), V = t(R) %*% R,log = FALSE){
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

mle.matrixnorm = function(data,row.restrict="none",col.restrict="none",tol = 1e-9,max.iter=100){
  # intend to implement toeplitz restriction later
  # if data is array, presumes indexed over third column
  # if list, presumes is a list of the matrices
  
}

