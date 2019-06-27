#   mixnorm.R
#   MixMatrix: Classification with Matrix Variate Normal and t distributions
#   Copyright (C) 2018-9  GZ Thompson <gzthompson@gmail.com>
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#   along with this program; if not, a copy is available at
#   https://www.R-project.org/Licenses/

##' Fit a matrix variate mixture model
##'
##' @param x data, \code{p x q x n} array
##' @param init a list containing an array of \code{K} of \code{p x q} means,
##'     and optionally \code{p x p} and \code{q x q} positive definite variance
##'     matrices. By default, those are presumed to be identity if not provided.
##'     If \code{init} is missing, it will be provided using the prior or K by
##'     \code{init_matrixmix}.
##' @param prior prior for the \code{K} classes, a vector that adds to unity
##' @param K number of classes - provide either this or the prior. If this is
##'     provided, the prior will be of equal distribution.
##' @param iter maximum number of iterations.
##' @param model whether to use the \code{normal} or \code{t} distribution.
##'     Currently, only the normal distribution is allowed.
##' @param method what method to use to fit the distribution. Currently no options.
##' @param tolerance convergence criterion, using Aitken acceleration of the
##'     log-likelihood.
##' @param nu degrees of freedom parameter
##' @param ... pass additional arguments to \code{MLmatrixnorm} or \code{MLmatrixt}
##' @param verbose whether to print diagnostic output, by default \code{0}. Higher
##'     numbers output more results.
##' @param miniter minimum number of iterations
##' @return A list of class \code{MixMatrixModel} containing the following
##'     components:
##' \describe{
##'      \item{\code{prior}}{the prior probabilities used.}
#'       \item{\code{K}}{the number of groups} 
#'       \item{\code{n}}{the number of observations} 
#'       \item{\code{centers}}{the group means.}
#'       \item{\code{U}}{the between-row covariance matrices}
#'       \item{\code{V}}{the between-column covariance matrix}
#'       \item{\code{posterior}}{the posterior probabilities for each observation}
#'       \item{\code{pi}}{ the final proportions}
#'       \item{\code{nu}}{The degrees of freedom parameter if the t distribution
#'            was used.}
#'       \item{\code{convergence }}{whether the model converged}
#'       \item{\code{logLik }}{a vector of the log-likelihoods of each iteration ending in
#'               the final log-likelihood of the model}
#'       \item{\code{method }}{the method used}
#'       \item{\code{call}}{The (matched) function call.}
##'    }
##'
##'
##' @export
##'
##' @examples
##'
##' set.seed(20180221)
#' A <- rmatrixnorm(30,mean=matrix(0,nrow=3,ncol=4))
#' B <- rmatrixnorm(30,mean=matrix(2,nrow=3,ncol=4))
#' C <- array(c(A,B), dim=c(3,4,60))
#' prior <- c(.5,.5)
#' init = list(centers = array(c(rep(0,12),rep(2,12)), dim = c(3,4,2)),
#'              U = array(c(diag(3), diag(3)), dim = c(3,3,2)),
#'              V = array(c(diag(4), diag(4)), dim = c(4,4,2))
#'              )
##'res<-matrixmixture(C, init = init, prior = prior)
##'print(res) # note: prints head of posterior, not full list
##'res<-matrixmixture(C, init = init, prior = prior, model = "t", nu = 5)
matrixmixture <- function(x, init = NULL, prior = NULL, K = length(prior), iter=1000,
                          model = "normal", method = NULL,
                          tolerance = 1e-1, nu=NULL, ..., verbose = 0, miniter = 3){
    if (class(x) == "list")
        x <- array(unlist(x),
                   dim = c(nrow(x[[1]]),
                           ncol(x[[1]]), length(x)))
    if (is.null(dim(x)))
        stop("'x' is not an array")
    if (any(!is.finite(x)))
        stop("infinite, NA or NaN values in 'x'")
    if (is.null(nu) || nu == 0 || is.infinite(nu)) model = "normal"
        
    if (model == "normal") nu = 0
    if (model != "normal") {
        df = nu
    }
        
    dims = dim(x)
    ## x is a p x q x n array
    n <- dims[3]
    p <- dims[1]
    q <- dims[2]
    if (verbose>0) cat("Dims: ",dims,"\n")
    if (!is.null(prior)) {
        if((length(prior) == 1) && (round(prior) == prior))
            prior = rep(1,prior)/prior
        
        if (any(prior < 0) || round(sum(prior), 5) != 1)
            stop("invalid 'prior'")
        prior <- prior[prior > 0L]
        K = length(prior)
    } else {
        if (missing(K)) stop("No prior and no K")

        prior = rep(1,K)/K
    }
    if(is.null(init))
        init = init_matrixmixture(x, prior = prior,...)

### extract initialization state
    ### should perhaps handle this by passing to init
    nclass = length(prior)
    centers = init$centers
    if( !is.null(init$U)){
        U = init$U
    } else {
        U = array(rep(diag(p),nclass),c(p,p,nclass))
    }
    if( !is.null(init$V)){
        V = init$V
    } else {
        V = array(rep(diag(q),nclass),c(q,q,nclass))
    }
    posterior = matrix(rep(prior, n),byrow = TRUE, nrow = n)
    newposterior = posterior
    eps = 1e40
    pi = prior
    logLikvec = numeric(0)
    if (verbose>1) {
        cat("\nInit centers: \n\n")
        print(init$centers)
        }
    if (verbose>2){
        print("Initial U and V")
        print(U)
        print(V)
    }
    convergeflag = FALSE
    Smatrix = array(0,c(p,p,n))
    SS = array(0,c(p,p,K))
    SSX = array(0,c(p,q,K))
    SSXX = array(0,c(q,q,K))
    newU = U
    newV = V
    newcenters = centers
    logLik = 0
    oldlogLik = 0
    olderlogLik = 0
    for(obs in 1:n){
        for(j in 1:K){
            logLik = logLik + log(pi[j]) +
                dmatrixt(x = x[,,obs], df = nu, mean = centers[,,j],
                         U = U[,,j], V = V[,,j], log = TRUE)
        }
    }
    logLikvec = c(logLikvec,logLik)
    
    i = 0
    while(i < iter && ( (abs(eps) > tolerance) || (i < miniter))){
        if(verbose) cat("\nEntering iteration:", i)
        if(verbose>1) print(pi)
        centers = newcenters
        newcenters = array(0, dim = c(p,q,K))
        U = newU
        V = newV
        posterior = newposterior
  
####### E STEP
        ## update expectations of sufficient statistics
        
        ## update z_ig weights
       # for(obs in 1:n){
            for(j in 1:K){
                newposterior[,j] = log(pi[j]) +
                    dmatrixt(x = x[,,],
                             df = nu, mean = centers[,,j],
                             U = U[,,j], V = V[,,j], log = TRUE)
              }
        #}
        newposterior <- ((newposterior - apply(newposterior, 1L, min, na.rm = TRUE)))
        newposterior = exp(newposterior)
        totalpost = rowSums(newposterior)
        newposterior = newposterior / totalpost
        if(verbose>1) print(newposterior[1:3,])
        
        ## update S_ig - conditional weights, only if non-normal
        
        if(model == "t"){
            dfmult = df + p + q - 1
            for(j in 1:K){
                ##### these don't work
                zigmult = rep(newposterior[,j], each = p*p)
                swept.data <- sweep(x, c(1, 2), centers[,,j])

           
                
                Stmp = xatx(swept.data,V[,,j])
                for (obs in 1:n) Stmp[,,obs] = Stmp[,,obs] + U[,,j]
                Smatrix = cubeinv(Stmp) * zigmult
                
                SS[,,j] = rowSums(Smatrix ,FALSE, 2)
                
                SSXtmp = cubemult(Smatrix, x)
                SSX[,,j] = rowSums(SSXtmp, FALSE, 2)
                
                SSXXtmp = cubemult(x,SSXtmp)
                SSXX[,,j] = rowSums(SSXXtmp,FALSE, 2)
                
            }
        }
        ### leave blank for now
        
####### CM STEPS
        pi = colMeans(newposterior)
        if (verbose) cat("\nNew pi: ", pi,"\n")
        ## max for centers, U, V
### max for centers
        ## if normal
        sumzig = colSums(newposterior)
        if(verbose>1) cat("\n Column sums of posterior", sumzig)
        if(model == "normal"){
            for(obs in 1:n){
                for(j in 1:K){
                    newcenters[,,j] = newcenters[,,j] + x[,,obs] * newposterior[obs,j]
                    if(verbose > 2) print(newcenters[,,j])
                }
            }
            for(j in 1:K) newcenters[,,j] = newcenters[,,j] / sumzig[j]
        } else {

            for(j in 1:K){
            newcenters[,,j] =  solve( SS[,,j]) %*% SSX[,,j]
            if(verbose > 2) print(newcenters[,,j])
            }
        }
                                        
        
    
### max for U, V
    ## if normal
    if(model == "normal"){
        for(j in 1:K){
            zigmult = rep(newposterior[,j], each = q*q)
            swept.data   <- sweep(x, c(1, 2), newcenters[,,j])
            inter.V <- txax(swept.data, U[,,j]) * zigmult
            newV[,,j] <- rowSums(inter.V, dims = 2)/(sumzig[j] * p)
            if(verbose >2) print(newV[,,j])
               
            zigmult = rep(newposterior[,j], each = p*p)
            inter.U <- xatx(swept.data, V[,,j]) * zigmult
            new.U = rowSums(inter.U, dims = 2)/(sumzig[j]*q)
            newU[,,j] <- new.U/(new.U[1, 1])
            if(verbose >2) print(newU[,,j])
        }
    } else {
        for(j in 1:K){

            newV[,,j] = (dfmult / (sumzig[j] * p)) * (SSXX[,,j] - t(SSX[,,j]) %*% newcenters[,,j] - t(newcenters[,,j]) %*% SSX[,,j] + t(newcenters[,,j]) %*% SS[,,j] %*% newcenters[,,j])
            newV[,,j] = newV[,,j]/newV[1,1,j]
            
            newUinv = (dfmult/(sumzig[j] * (df + p - 1))) * SS[,,j]
            newU[,,j] = solve(newUinv)
        }
    }
        
####### Eval convergence
        if(verbose > 1){
            print("New centers:")
            print(newcenters)
            print("New U:")
            print(U)
            print("New V:")
            print(V)
            }

        olderlogLik = oldlogLik
        oldlogLik = logLik
        logLik = 0
        for(obs in 1:n){
            for(j in 1:K){
                logLik = logLik + log(pi[j]) +
                dmatrixt(x = x[,,obs], df = nu, mean = newcenters[,,j],
                         U = newU[,,j], V = newV[,,j], log = TRUE)
            }
        }
        if(verbose) cat("\nLog likelihood:", logLik)
        if(i == 0) {
            ## initialize to some not-so-bad values so that doesn't immediately "converge"
            olderlogLik = oldlogLik - .2*abs(oldlogLik)
            }
        #eps = sum((newcenters - centers)^2)+sum( (newU-U)^2) + sum( (newV-V)^2 )
        aitken = (logLik - oldlogLik) / (oldlogLik - olderlogLik)
        linf = oldlogLik + 1/(1-aitken) * (logLik - oldlogLik)
        eps = linf - logLik
        i = i + 1
        if(verbose) cat("\nAitken, l_infinity, epsilon:", aitken, linf, eps)
        logLikvec = c(logLikvec, logLik)
    }
    if ((i == iter || eps > tolerance) ){
        warning('failed to converge')
    } else convergeflag <- TRUE
    if(verbose) cat("\nDone at iteration ", i-1,"\n")
    U = newU
    V = newV
    centers = newcenters
    posterior = newposterior
    pi = colMeans(posterior)
    logLik = 0
    for(i in 1:n){
        for(j in 1:K){
            logLik = logLik + log(pi[j]) +
                dmatrixt(x = x[,,i], df = nu, mean = centers[,,j],
                         U = U[,,j], V = V[,,j], log = TRUE)
        }
    }
    if(verbose>1) {
        print("Final centers:")
        print(centers)
        }
    if(verbose) cat("\nLog Likelihood Trace: \n", logLikvec, "\n")
    cl <- match.call()
    cl[[1L]] <- as.name("matrixmixture")

    structure(
    list(
        prior = prior,
        K = nclass,
        n = n,
        centers = centers,
        U = U,
        V = V,
        posterior = posterior,
        pi = pi,
        convergence = convergeflag,
        iter = i, 
        logLik = logLikvec,
        model = model,
        method = method,
        call = cl
    ),
    class = "MixMatrixModel")
}


print.MixMatrixModel <- function(model){
    model[["posterior"]] = head(model[["posterior"]])             
    print.default(model)
}

##' Initializing settings for Matrix Mixture Models
##'
##' Providing this will generate a list suitable for use as the \code{init}
##' argument in the \code{matrixmixture} function. Either provide data
##' and it will select centers and variance matrices to initialize or
##' provide initial values and it will format them as expected for the function.
##' 
##' @param data array of data
##' @param prior prior probability. One of \code{prior} and \code{K} must be
##'      provided. They must be consistent if both provided.
##' @param K number of groups
##' @param centers either a matrix or an array of matrices for use as the
##'      \code{centers} argument (optional)
##' @param U either a matrix or an array of matrices for use as the \code{U}
##'      argument (optional)
##' @param V either a matrix or an array of matrices for use as the \code{V}
##'      argument (optional)
##' @param centermethod what method to use to generate initial centers.
##'      Current support random start or performing $k$-means on the vectorized
##'      version for a small number of iterations and then converting back.
##'      By default, if centers are provided, nothing will be done.
##' @param varmethod what method to use to choose initial variance matrices.
##'      Currently either identity matrices or the empirical covariance matrix
##'      determined by hard assignment to the nearest centers.
##'      By default, if \code{U} and \code{V} matrices are provided, nothing
##'      will be done.
##' @param model whether to use a normal distribution or a t-distribution, not
##'      relevant for more initialization methods.
##' @param ... Additional arguments to pass to $k$-means.
##' @return a list suitable to use as the \code{init} argument in
##'      \code{matrixmixture}
##' @importFrom stats kmeans
init_matrixmixture<- function(data, prior, K = length(prior), centers = NULL,
                              U = NULL, V = NULL,  centermethod = "kmeans",
                              varmethod = "identity", model = "normal",...){
    dims = dim(data)
    p = dims[1]
    q = dims[2]
    n = dims[3]

    if(length(prior) == 1) K = prior
    
    if(centermethod == "random"){
    select = sample(n,K, replace = FALSE)
    centers = data[,,select]
    }
    if(centermethod == "kmeans" || centermethod == "k-means"){
        res = kmeans(matrix(data, nrow = n), centers = K, ...)
        centers = array(res$centers, dim = c(p,q,K))
    }
    if(!is.null(U) && !is.null(V)){
        
        if(length(dim(U) == 2)) U = array(rep(U,K), dim = c(p,p,K))
        if(length(dim(V) == 2)) V = array(rep(V,K), dim = c(q,q,K))
        
    } else {
        if(varmethod == "identity"){
            U = array(c(rep(diag(p),K)), dim = c(p,p,K))
            V = array(c(rep(diag(q),K)), dim = c(q,q,K))
        }
    }
    list(
        centers = centers,
        U = U,
        V = V
        )
}
