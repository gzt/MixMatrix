rmatrixnorm <-
function(n, mean,L=diag(dim(mean)[1]),R=diag(dim(mean)[2]),
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
