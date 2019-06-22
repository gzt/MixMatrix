#   mixnorm.R
#   MixMatrix: Classification with Matrix Variate Normal and t distributions
#   Copyright (C) 2018-9  GZ Thompson <gzthompson@gmail.com>
#
#   These functions are based on modifications of the source for 
#   MASS::lda() and MASS::qda(),
#   copyright (C) 1994-2013 W. N. Venables and B. D. Ripley
#   released under GPL 2 or greater. This software is released under GPL-3.
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
##' @title 
##' @param x data, \code{p x q x n} array
##' @param prior prior for the \code{K} classes, a vector that adds to unity
##' @param init a list containing an array of \code{K} of \code{p x q} means,
##'     and optionally \code{p x p} and \code{q x q} positive definite variance
##'     matrices. By default, those are presumed to be identity if not provided.
##' @param iter maximum number of iterations.
##' @param method whether to use the \code{normal} or \code{t} distribution.
##'     Currently, only the normal distribution is allowed.
##' @param tol convergence criterion
##' @param nu degrees of freedom parameter
##' @param ... pass additional arguments to \code{MLmatrixnorm} or \code{MLmatrixt}
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
#'       \item{\code{logLik }}{the final log-likelihood of the model}
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
#' B <- rmatrixnorm(30,mean=matrix(1,nrow=3,ncol=4))
#' C <- array(c(A,B), dim=c(3,4,60))
#' prior <- c(.5,.5)
#' init = list(centers = array(c(rep(0,12),rep(1,12)), dim = c(3,4,2)),
#'              U = array(c(diag(3), diag(3)), dim = c(3,3,2)),
#'              V = array(c(diag(4), diag(4)), dim = c(4,4,2))
#'              )
##'matrixmixture(C, prior, init)
##'
##' 
matrixmixture <- function(x, prior, init, iter=1000, method = "normal", tol = 1e-6, nu=NULL, ...){

    if (class(x) == "list")
        x <- array(unlist(x),
                   dim = c(nrow(x[[1]]),
                           ncol(x[[1]]), length(x)))
    if (is.null(dim(x)))
        stop("'x' is not an array")
if (any(!is.finite(x)))
    stop("infinite, NA or NaN values in 'x'")
    if (nu == 0 || is.infinite(nu)) method = "normal"
    
if (method == "normal") nu = NULL
    if (method != "normal") {
        method = "normal"
        warning("t not implemented yet, using normal")
    }

    dims = dim(x)
    # x is a p x q x n array
    n <- dims[3]
    p <- dims[1]
    q <- dims[2]
    if (!missing(prior)) {
        if (any(prior < 0) || round(sum(prior), 5) != 1)
            stop("invalid 'prior'")
        prior <- prior[prior > 0L]
    }
    
### take in data
    convergeflag = FALSE
    dims = dim(data)
### extract initialization state
    nclass = length(prior)
    centers = init$centers
    if( !is.null(init$U)){
        U = init$U
    } else {
        U = array(rep(diag(p),nclass),c(p,p,nclass))
    }
    if( !is.null(init$V)){
        U = init$V
    } else {
        V = array(rep(diag(q),nclass),c(q,q,nclass))
    }
    posterior = matrix(rep(prior, n),byrow = TRUE, nrow = n)
    eps = 1e40
    
    while(i < iter && eps > tol){
        newcenters = centers
        newU = U
        newV = V
        newposterior = posterior
####### E STEP
        ## update expectations of sufficient statistics
        
        ## update pi_i weights

####### M STEP
        ## max for centers, U, V


####### Eval convergence
        eps = sum((newcenters - centers)^2)+sum( (newU-U)^2) + sum( (newV-V)^2 )
        U = newU
        V = newV
        centers = newcenters
        posterior = newposterior
    }
    if (i == iter || eps > tol){
        warning('failed to converge')
        }


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
        pi = colMeans(posterior),
        convergence = convergeflag,
        logLik = logLik,
        method = method,
        call = cl
    ),
    class = "MixMatrixModel")
}
