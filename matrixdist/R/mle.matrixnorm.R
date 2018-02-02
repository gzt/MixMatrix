mle.matrixnorm <-
function(data,row.restrict="none",col.restrict="none",tol = 1e-9,max.iter=100,U,V){
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
    error.term = max(sum((new.V - V)^2), sum((new.U - U)^2))
    V = new.V
    U = new.U
    iter = iter + 1
  }
  if(iter >= max.iter || error.term > tol) warning("Failed to converge")
  
  return(list(mean = mu, U=U, V = V, iter=iter, tol=error.term))
}
