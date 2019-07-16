

.SStep <- function(data,centers,U,V,weights){
    dims <- dim(data)
    p <- dims[1]
    q <- dims[2]
    n <- dims[3]
    
    
    zigmult = rep(weights, each = p*p)
    swept.data <- sweep(data, c(1, 2), centers)
    
    Stmp = xatx(swept.data,V)
    for (obs in 1:n) Stmp[,,obs] = Stmp[,,obs] + U
    Smatrix = cubeinv(Stmp) * zigmult
    
    SS = rowSums(Smatrix ,FALSE, 2)
    
    SSXtmp = cubemult(Smatrix, data)
    SSX = rowSums(SSXtmp, FALSE, 2)
    
    SSXXtmp = cubemult(data,SSXtmp)
    SSXX = rowSums(SSXXtmp,FALSE, 2)
    SSD = detsum(Smatrix)

    list(SS=SS, SSX = SSX, SSXX = SSXX, SSD = SSD)
}

.MeansFunction <- function(data, U=NULL,V=NULL, SS, SSX, weights, row.mean, col.mean, model){
    dims <- dim(data)
    p <- dims[1]
    q <- dims[2]
    n <- dims[3]
    
    sumzig = sum(weights)
    newcenters = matrix(0,nrow=p,ncol=q)
    if(model == "normal"){
        for(obs in 1:n){
            newcenters = newcenters + data[,,obs] * weights[obs]
        }
        
        newcenters = newcenters / sumzig
        if (row.mean) {
                                        # make it so that the mean is constant within a row
            newcenters <- matrix(rowMeans(newcenters), nrow = dims[1], ncol = dims[2])
        }
        if (col.mean) {
                                        # make it so that the mean is constant within a column
            newcenters <- matrix(colMeans(newcenters), nrow = dims[1], ncol = dims[2], byrow = TRUE)
        }
        
    } else {
        
        
        if (row.mean && col.mean) {
                                        # make it so that the mean is constant within a row
            scalarmu = matrixtrace(SSX %*% solve(V) %*% ones(q,p)) / matrixtrace(SS %*% ones(p,q) %*% solve(V) %*% ones(q,p))
            newcenters <-   scalarmu * ones(p,q)
        } else if (col.mean) {
                                        # make it so that the mean is constant within a column
                                        # ie mu = p x 1, times ones 1 x q
            newcenters <- ones(p,p) %*% SSX / sum(SS)
        } else if (row.mean) {
                                        # make it so that the mean is constant within a row
                                        # ie  ones p x 1 times mu = 1 x q
            newcenters = solve( SS) %*% SSX %*% (solve(V) %*% ones(q,q)) / sum(solve(V))
        } else {
            newcenters =  solve( SS) %*% SSX
        }
    }
    
    newcenters
}





.colVars <- function(data,U,V,df,SS,SSX,SSXX,col.variance,col.set.var){

}


.rowVars <- function(data,U,V,df,SS,SSX,SSXX,col.variance,col.set.var){

}


.varparse <- function(varoption){
    varflag = FALSE
    varopt = varoption
    
    if (grepl("^i", x = varoption,ignore.case = TRUE)) {
        varflag = TRUE
        varopt = "I"
    }
    
    if (grepl("^co", x = varoption,ignore.case = TRUE)) {
        varflag = FALSE
        varopt = "cor"
    }
    if (grepl("^ar", x = varoption,ignore.case = TRUE)) {
        varflag = TRUE 
        varopt= "AR(1)"
    }
    if (grepl("^cs", x = varoption,ignore.case = TRUE)) {
        varflag = TRUE
        varopt = "CS"
    }
    
    list(varflag=varflag,varopt=varopt)
}
