dmatrixnorm <-
function(x, mean = array(0L, dim(x)), L = diag(dim(x)[1]), R=diag(dim(x)[2]), 
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
