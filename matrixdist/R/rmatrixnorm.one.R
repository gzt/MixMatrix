rmatrixnorm.one <-
function(mean,L=diag(dim(mean)[1]),R=diag(dim(mean)[2]), U = L %*% t(L), V = t(R) %*% R){
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
