dmatrixnorm.test <-
function(x, mean = array(0L, dim(x)), L = diag(dim(x)[1]), R=diag(dim(x)[2]), 
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
