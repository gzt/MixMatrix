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
##' Clustering by fitting a mixture model using EM with \code{K} groups
##' and unconstrained covariance matrices for a matrix variate normal or
##' matrix variate t distribution (with specified degrees of freedom \code{nu}). 
##'
##' @param x data, \eqn{p \times q \times n}{p * q * n} array
##' @param init a list containing an array of \code{K} of \eqn{p \times q}{p * q} means labeled
##'     \code{centers},
##'     and optionally \eqn{p \times p}{p * p} and \eqn{q \times q}{q * q} positive definite variance
##'     matrices labeled \code{U} and \code{V}.
##'     By default, those are presumed to be identity if not provided.
##'     If \code{init} is missing, it will be provided using the \code{prior} or \code{K} by
##'     \code{init_matrixmix}.
##' @param prior prior for the \code{K} classes, a vector that adds to unity
##' @param K number of classes - provide either this or the prior. If this is
##'     provided, the prior will be of equal distribution.
##' @param iter maximum number of iterations.
##' @param model whether to use the \code{normal} or \code{t} distribution.
##'     
##' @param method what method to use to fit the distribution. Currently no options.
##' @param row.mean By default, \code{FALSE}. If \code{TRUE}, will fit a
##'    common mean within each row. If both this and \code{col.mean} are
##'    \code{TRUE}, there will be a common mean for the entire matrix.
##' @param col.mean By default, \code{FALSE}. If \code{TRUE}, will fit a
##'    common mean within each row. If both this and \code{row.mean} are
##'    \code{TRUE}, there will be a common mean for the entire matrix.
##' @param tolerance convergence criterion, using Aitken acceleration of the
##'     log-likelihood by default.
##' @param nu degrees of freedom parameter
##' @param ... pass additional arguments to \code{MLmatrixnorm} or \code{MLmatrixt}
##' @param verbose whether to print diagnostic output, by default \code{0}. Higher
##'     numbers output more results.
##' @param miniter minimum number of iterations
##' @param convergence By default, \code{TRUE}. Whether to use Aitken acceleration to
##'     determine convergence. If false, it instead checks if the change in
##'     log-likelihood is less than \code{tolerance}.
##' @return A list of class \code{MixMatrixModel} containing the following
##'     components:
##' \describe{
##'      \item{\code{prior}}{the prior probabilities used.}
##'      \item{\code{init}}{the initialization used.}
#'       \item{\code{K}}{the number of groups} 
#'       \item{\code{N}}{the number of observations} 
#'       \item{\code{centers}}{the group means.}
#'       \item{\code{U}}{the between-row covariance matrices}
#'       \item{\code{V}}{the between-column covariance matrix}
#'       \item{\code{posterior}}{the posterior probabilities for each observation}
#'       \item{\code{pi}}{ the final proportions}
#'       \item{\code{nu}}{The degrees of freedom parameter if the t distribution
#'            was used.}
#'       \item{\code{convergence }}{whether the model converged}
#'       \item{\code{logLik}}{a vector of the log-likelihoods
#'               of each iteration ending in
#'               the final log-likelihood of the model}
#'       \item{\code{model}}{the model used}
#'       \item{\code{method}}{the method used}
#'       \item{\code{call}}{The (matched) function call.}
##'    }
##'
##'
##' @export
##' @seealso \code{\link{init_matrixmixture}}
##'
##' @references
##'     Andrews, Jeffrey L., Paul D. McNicholas, and Sanjeena Subedi. 2011.
##'       "Model-Based Classification via Mixtures of Multivariate
##'       T-Distributions." Computational Statistics & Data Analysis 55 (1):
##'       520–29. \doi{10.1016/j.csda.2010.05.019}.
##'
##'     Fraley, Chris, and Adrian E Raftery. 2002. "Model-Based Clustering,
##'        Discriminant Analysis, and Density Estimation." Journal of the
##'        American Statistical Association 97 (458). Taylor & Francis: 611–31.
##'        \doi{10.1198/016214502760047131}.
##'
##'     McLachlan, Geoffrey J, Sharon X Lee, and Suren I Rathnayake. 2019.
##'           "Finite Mixture Models." Annual Review of Statistics and Its
##'           Application 6. Annual Reviews: 355–78.
##'           \doi{10.1146/annurev-statistics-031017-100325}.
##'
##'     Viroli, Cinzia. 2011. "Finite Mixtures of Matrix Normal Distributions
##'           for Classifying Three-Way Data." Statistics and Computing 21 (4):
##'           511–22. \doi{10.1007/s11222-010-9188-x}.
##' 
##' @examples
##' set.seed(20180221)
##' A <- rmatrixt(30,mean=matrix(0,nrow=3,ncol=4), df = 5)
##' # 3x4 matrices with mean 0
##' B <- rmatrixt(30,mean=matrix(1,nrow=3,ncol=4), df = 5)
##' # 3x4 matrices with mean 2
##' C <- array(c(A,B), dim=c(3,4,60)) # combine into one array
##' prior <- c(.5,.5) # equal probability prior
##' # create an intialization object, starts at the true parameters
##' init = list(centers = array(c(rep(0,12),rep(1,12)), dim = c(3,4,2)),
##'               U = array(c(diag(3), diag(3)), dim = c(3,3,2))*20,
##'               V = array(c(diag(4), diag(4)), dim = c(4,4,2))
##'  )
##' # fit model
##'  res<-matrixmixture(C, init = init, prior = prior, nu = 5,
##'                     model = "t", tolerance = 1e-1)
##' print(res$centers) # the final centers
##' print(res$pi) # the final mixing proportion
##' plot(res) # the log likelihood by iteration
##' logLik(res)
##' BIC(res)
##' predict(res, newdata = C[,,c(1,31)])
matrixmixture <- function(x, init = NULL, prior = NULL, K = length(prior), iter=1000,
                          model = "normal", method = NULL, row.mean = FALSE, col.mean = FALSE,
                          tolerance = 1e-1, nu=NULL, ..., verbose = 0, miniter = 5, convergence = TRUE){
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
    ## x is a p * q * n array
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
        if(model == "t") U = (nu-2) * stats::var(x[1,1,]) * U
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
    i = 0
    while(i < iter && ( ((eps) > tolerance) || (i < miniter))){
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

            for(j in 1:K){
                newposterior[,j] = log(pi[j]) +
                    dmatrixt(x = x[,,],
                             df = nu, mean = centers[,,j],
                             U = U[,,j], V = V[,,j], log = TRUE)
              }

        newposterior <- ((newposterior - apply(newposterior, 1L, min, na.rm = TRUE)))
        newposterior = exp(newposterior)
        totalpost = rowSums(newposterior)
        newposterior = newposterior / totalpost
        if(verbose>1) print(newposterior[1:5,])
        
        ## update S_ig - conditional weights, only if non-normal
        
        if(model == "t"){
            dfmult = df + p + q - 1
            for(j in 1:K){
                
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
            for(j in 1:K) {
                newcenters[,,j] = newcenters[,,j] / sumzig[j]
                if (row.mean) {
                                        # make it so that the mean is constant within a row
                    newcenters[,,j] <- matrix(rowMeans(newcenters[,,j]), nrow = dims[1], ncol = dims[2])
                }
                if (col.mean) {
                                        # make it so that the mean is constant within a column
                    newcenters[,,j] <- matrix(colMeans(newcenters[,,j]), nrow = dims[1], ncol = dims[2], byrow = TRUE)
                }
            }
        } else {
            
            for(j in 1:K){
                if (row.mean && col.mean) {
                                        # make it so that the mean is constant within a row
                    scalarmu = matrixtrace(SSX[,,j] %*% solve(V[,,j]) %*% ones(q,p)) / matrixtrace(SS[,,j] %*% ones(p,q) %*% solve(V[,,j]) %*% ones(q,p))
                    newcenters[,,j] <-   scalarmu * ones(p,q)
                } else if (col.mean) {
                                        # make it so that the mean is constant within a column
                                        # ie mu = p x 1, times ones 1 x q
                    newcenters[,,j] <- ones(p,p) %*% SSX[,,j] / sum(SS[,,j])
                } else if (row.mean) {
                                        # make it so that the mean is constant within a row
                                        # ie  ones p x 1 times mu = 1 x q
                    newcenters[,,j] = solve( SS[,,j]) %*% SSX[,,j] %*% (solve(V[,,j]) %*% ones(q,q)) / sum(solve(V[,,j]))
                } else {
                    newcenters[,,j] =  solve( SS[,,j]) %*% SSX[,,j]
                }
                if(verbose > 2) print(newcenters[,,j])
            }
        }
            
### max for U, V
    ## if normal
    if(model == "normal"){
        for(j in 1:K){
           # for (iterations in 1:20){
            zigmult = rep(newposterior[,j], each = q*q)
            swept.data   <- sweep(x, c(1, 2), newcenters[,,j])
            inter.V <- txax(swept.data, U[,,j]) * zigmult
            newV[,,j] <- rowSums(inter.V, dims = 2)/(sumzig[j] * p)
            if(verbose >2) print(newV[,,j])
               
            zigmult = rep(newposterior[,j], each = p*p)
            inter.U <- xatx(swept.data, newV[,,j]) * zigmult
            new.U = rowSums(inter.U, dims = 2)/(sumzig[j]*q)
            newU[,,j] <- new.U/(new.U[1, 1])
            if(verbose >2) print(newU[,,j])
           # }
        }
    } else {
        for(j in 1:K){

            newV[,,j] = (dfmult / (sumzig[j] * p)) * (SSXX[,,j] - t(SSX[,,j]) %*% newcenters[,,j] -
                                                      t(newcenters[,,j]) %*% SSX[,,j] + t(newcenters[,,j]) %*% SS[,,j] %*% newcenters[,,j])
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
            print(newU)
            print("New V:")
            print(newV)
            }

        olderlogLik = oldlogLik
        oldlogLik = logLik
        logLik = 0
        for(obs in 1:n){
            for(j in 1:K){
                logLik = logLik + newposterior[obs,j]*(log(pi[j]) +
                dmatrixt(x = x[,,obs], df = nu, mean = newcenters[,,j],
                         U = newU[,,j], V = newV[,,j], log = TRUE))
            }
        }
        if(verbose) cat("\nLog likelihood:", logLik)
        if(i == 0) {
            ## initialize to some not-so-bad values so that doesn't immediately "converge"
            olderlogLik = oldlogLik - .2*abs(oldlogLik)
            }
                                        #eps = sum((newcenters - centers)^2)+sum( (newU-U)^2) + sum( (newV-V)^2 )
        if(convergence){
        aitken = (logLik - oldlogLik) / (oldlogLik - olderlogLik)
        linf = oldlogLik + 1/(1-aitken) * (logLik - oldlogLik)
        eps = linf - logLik
        } else eps = logLik - oldlogLik
        
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
        init = init,
        K = nclass,
        N = n,
        centers = centers,
        U = U,
        V = V,
        posterior = posterior,
        pi = pi,
        nu = nu,
        convergence = convergeflag,
        iter = i, 
        logLik = logLikvec,
        model = model,
        method = method,
        call = cl
    ),
    class = "MixMatrixModel")
}

#' @export
print.MixMatrixModel <- function(x, ...){
    x[["posterior"]] = head(x[["posterior"]])
    x[["init"]] = NULL
    print.default(x,...)
    
}

#' @export
#' @importFrom graphics plot
plot.MixMatrixModel <- function(x, ...){
    plot(x = 1:length(x$logLik), y = x$logLik,
         ylab="Log Likelihood",xlab="iteration",...)
}


#' @export
logLik.MixMatrixModel <- function(object, ...){

    dims = dim(object$centers)
    n <- object$call$N
    p <- dims[1]
    q <- dims[2]
    numgroups  = length(levels(grouping))
    grouplist = levels(grouping)
    meanpars = p*q
    if(!is.null(object$call$row.mean) && (object$call$row.mean)) meanpars = meanpars/q
    if(!is.null(object$call$col.mean) && (object$call$col.mean)) meanpars = meanpars/p
    upars = (p+1)*p/2
    vpars = (q+1)*q/2 # note of course that there's one par that will get subbed off variance
    nupar = 0 # only fixed for now
    numgroups = (object$K)

### insert here logic for parsing out different values for this later
### as ways of restricting variances and means are added
    
    df = numgroups*(vpars + upars + nupar + meanpars - 1)
    logLik = object$logLik[length(object$logLik)]
    
    class(logLik) = "logLik"
    attr(logLik, 'df') <- df
    attr(logLik, 'nobs') <- n
    logLik
    
}

#' @export
nobs.MixMatrixModel <- function(object, ...){
    object$N
}

##' Initializing settings for Matrix Mixture Models
##'
##' Providing this will generate a list suitable for use as the \code{init}
##' argument in the \code{matrixmixture} function. Either provide data
##' and it will select centers and variance matrices to initialize or
##' provide initial values and it will format them as expected for the function.
##' 
##' @param data  data, \eqn{p \times q \times n}{p * q * n} array
##' @param prior prior probability. One of \code{prior} and \code{K} must be
##'      provided. They must be consistent if both provided.
##' @param K number of groups
##' @param centers (optional) either a matrix or an array of \eqn{p \times p}{p * p}
##'      matrices for use as the \code{centers} argument.
##'      If fewer than \code{K} are provided, the
##'      remainder are chosen by \code{centermethod}.
##' @param U (optional) either a matrix or an array of  \eqn{q \times q}{q * q}
##'      matrices for use as the \code{U}
##'      argument. If a matrix is provided, it is duplicated to provide an array.
##'      If an array is provided, it should have \code{K} slices.
##' @param V  (optional) either a matrix or an array of matrices for use as the \code{U}
##'      argument. If a matrix is provided, it is duplicated to provide an array.
##'      If an array is provided, it should have \code{K} slices.
##' @param centermethod what method to use to generate initial centers.
##'      Currently support random start (\code{random}) or performing k-means
##'      (\code{kmeans}) on the vectorized version for a small number of
##'      iterations and then converting back.
##'      By default, if \code{K} centers are provided, nothing will be done.
##' @param varmethod what method to use to choose initial variance matrices.
##'      Currently only identity matrices are created. 
##'      By default, if \code{U} and \code{V} matrices are provided, nothing
##'      will be done.
##' @param model whether to use a normal distribution or a t-distribution, not
##'      relevant for more initialization methods.
##' @param init (optional) a (possibly partially-formed) list with some of the components
##'     \code{centers}, \code{U}, and \code{V}. The function will complete the
##'     list and fill out missing entries.
##' @param ... Additional arguments to pass to \code{kmeans()} if that is
##'     \code{centermethod}.
##' @return a list suitable to use as the \code{init} argument in
##'      \code{matrixmixture}:
##' \describe{
#'       \item{\code{centers}}{the group means,
#'             a \eqn{p \times q \times K}{p * q * K} array.}
#'       \item{\code{U}}{the between-row covariance matrices, a \eqn{p \times p \times K}{p * p * K}  array}
#'       \item{\code{V}}{the between-column covariance matrix, a \eqn{q \times q \times K}{q * q * K} array}
##'    }
##'
##' @export
##' @importFrom stats kmeans
##' @seealso \code{\link{matrixmixture}}
##' 
##' @examples 
##'  set.seed(20180221)
##' A <- rmatrixt(30,mean=matrix(0,nrow=3,ncol=4), df = 10)
##' # 3x4 matrices with mean 0
##' B <- rmatrixt(30,mean=matrix(2,nrow=3,ncol=4), df = 10)
##' # 3x4 matrices with mean 2
##' C <- array(c(A,B), dim=c(3,4,60)) # combine into one array
##' prior <- c(.5,.5) # equal probability prior
##' init = init_matrixmixture(C, prior = prior)
##' # will find two centers using the "kmeans" method on the vectorized matrices
init_matrixmixture<- function(data, prior = NULL, K = length(prior), centers = NULL,
                              U = NULL, V = NULL,  centermethod = "kmeans",
                              varmethod = "identity", model = "normal", init = NULL,...){
    dims = dim(data)
    p = dims[1]
    q = dims[2]
    n = dims[3]
    remains = K
    cenflag = FALSE
    centers = array(dim = c(p,q,K))
    if(length(prior) == 1) K = prior
    if(!is.null(K)) prior = rep(1,K)/K
    if(!is.null(init)){
        if(!is.null(init$centers)){
            cenflag = TRUE
            initcenters = init$centers
            dimcen = dim(initcenters)
            if(!((dimcen[1]==p)&&(dimcen[2] == q))) stop("wrong dimension for provided centers")
            if(length(dimcen) == 2) remains = K -1
                else remains = K - dimcen[3]
        }
        if(is.null(U)) U = init$U
        if(is.null(V)) V = init$V
    }
    if(cenflag) centers[,,(remains+1):K] = initcenters
    
    if(centermethod == "random" && (remains > 0)){
        select = sample(n,remains, replace = FALSE)
        centers[,,1:remains] = data[,,select]
    }
    if((remains > 0) && (centermethod == "kmeans" || centermethod == "k-means")){
        res = kmeans(matrix(data, nrow = n), centers = remains, ...)
        centers = array(res$centers, dim = c(p,q,remains))
        if(cenflag) centers = array(c(centers, initcenters), dim = c(p,q,K))
    }
    if(!is.null(U)){
        if(length(dim(U) == 2)) U = array(rep(U,K), dim = c(p,p,K))
    }
    if(!is.null(V)){        
        if(length(dim(V) == 2)) V = array(rep(V,K), dim = c(q,q,K))
    } 
    if(varmethod == "identity"){
        if(is.null(U)) U = array(c(rep(diag(p),K)), dim = c(p,p,K)) 
        if(is.null(V)) V = array(c(rep(diag(q),K)), dim = c(q,q,K))
    }
   
    list(
        centers = centers,
        U = U,
        V = V
        )
}

#' @export
# S3 method for predict on class MixMatrixModel
predict.MixMatrixModel <- function(object, newdata, prior = object$prior,...){
        if (!inherits(object, "MixMatrixModel"))
            stop("object not of class \"MixMatrixModel\"")

    if (missing(newdata)) {
        newdata <- eval.parent(object$call$x)
    }

    if (is.null(dim(newdata)))
      stop("'newdata' is not an array")
    if (any(!is.finite(newdata)))
      stop("infinite, NA or NaN values in 'newdata'")
    x <- (newdata)
  

    if (length(dim(x)) == 2) x <- array(x, dim= c(dim(x),1))

    if (ncol(x[, , 1, drop = FALSE]) != ncol(object$centers[, , 1, drop = FALSE]))
      stop("wrong column dimension of matrices")
    if (nrow(x[, , 1, drop = FALSE]) != nrow(object$centers[, , 1, drop = FALSE]))
      stop("wrong row dimension of matrices")
    ng <- length(object$prior)
    if (!missing(prior)) {
      if (length(prior) != ng) stop("invalid prior length")
      if (any(prior < 0) || round(sum(prior), 5) != 1)
        stop("invalid 'prior'")
    }

    dims = dim(x)
    # x is a p x q x n array
    n <- dims[3]
    p <- dims[1]
    q <- dims[2]
    df <- object$nu

 
  dist = matrix(0, nrow = n, ncol = ng)
    posterior = matrix(0, nrow = n, ncol = ng)
   

    for (i in seq(n)) {
      Xi = matrix(x[, , i], p, q)
      for (j in seq(ng)) {
          dist[i,j] = log(prior[j]) + dmatrixt(x = Xi, df = df, mean = matrix(object$centers[,,j],nrow=p,ncol=q),
                               U = matrix(object$U[,,j],nrow=p,ncol=p), V = matrix(object$V[,,j],nrow=q,ncol=q), log = TRUE)
      }
    }
    posterior = exp( (dist - apply(dist, 1L, max, na.rm = TRUE)))
    totalpost = rowSums(posterior)
    posterior = posterior / totalpost
    nm <- names(object$prior)
    cl <- max.col(posterior)
    list(class = cl, posterior = posterior)      
        
}
