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
##' @param init a list containing an array of \code{K} of
##'     \eqn{p \times q}{p * q} means labeled \code{centers},
##'     and optionally \eqn{p \times p}{p * p} and \eqn{q \times q}{q * q}
##'     positive definite variance matrices labeled \code{U} and \code{V}.
##'     By default, those are presumed to be identity if not provided.
##'     If \code{init} is missing, it will be provided using the \code{prior}
##'     or \code{K} by \code{init_matrixmix}.
##' @param prior prior for the \code{K} classes, a vector that adds to unity
##' @param K number of classes - provide either this or the prior. If this is
##'     provided, the prior will be of uniform distribution among the classes.
##' @param iter maximum number of iterations.
##' @param model whether to use the \code{normal} or \code{t} distribution.
##'
##' @param method what method to use to fit the distribution.
##'     Currently no options.
##' @param row.mean By default, \code{FALSE}. If \code{TRUE}, will fit a
##'    common mean within each row. If both this and \code{col.mean} are
##'    \code{TRUE}, there will be a common mean for the entire matrix.
##' @param col.mean By default, \code{FALSE}. If \code{TRUE}, will fit a
##'    common mean within each row. If both this and \code{row.mean} are
##'    \code{TRUE}, there will be a common mean for the entire matrix.
##' @param tolerance convergence criterion, using Aitken acceleration of the
##'     log-likelihood by default.
##' @param nu degrees of freedom parameter. Can be a vector of length \code{K}.
##' @param ... pass additional arguments to \code{MLmatrixnorm} or
##'     \code{MLmatrixt}
##' @param verbose whether to print diagnostic output, by default \code{0}.
##'     Higher  numbers output more results.
##' @param miniter minimum number of iterations
##' @param convergence By default, \code{TRUE}, using Aitken acceleration
##'     to determine convergence. If false, it instead checks if the change in
##'     log-likelihood is less than \code{tolerance}. Aitken acceleration may
##'     prematurely end in the first few steps, so you may wish to set
##'     \code{miniter} or select \code{FALSE} if this is an issue.
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
#'       \item{\code{posterior}}{the posterior probabilities for each
#'             observation}
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
##' A <- rmatrixt(20,mean=matrix(0,nrow=3,ncol=4), df = 5)
##' # 3x4 matrices with mean 0
##' B <- rmatrixt(20,mean=matrix(1,nrow=3,ncol=4), df = 5)
##' # 3x4 matrices with mean 1
##' C <- array(c(A,B), dim=c(3,4,40)) # combine into one array
##' prior <- c(.5,.5) # equal probability prior
##' # create an intialization object, starts at the true parameters
##' init = list(centers = array(c(rep(0,12),rep(1,12)), dim = c(3,4,2)),
##'               U = array(c(diag(3), diag(3)), dim = c(3,3,2))*20,
##'               V = array(c(diag(4), diag(4)), dim = c(4,4,2))
##'  )
##' # fit model
##'  res<-matrixmixture(C, init = init, prior = prior, nu = 5,
##'                     model = "t", tolerance = 1e-3, convergence = FALSE)
##' print(res$centers) # the final centers
##' print(res$pi) # the final mixing proportion
##' plot(res) # the log likelihood by iteration
##' logLik(res) # log likelihood of final result
##' BIC(res) # BIC of final result
##' predict(res, newdata = C[,,c(1,21)]) # predicted class membership
matrixmixture <- function(x, init = NULL, prior = NULL, K = length(prior),
                          iter = 1000, model = "normal", method = NULL,
                          row.mean = FALSE, col.mean = FALSE,
                          tolerance = 1e-1, nu = NULL, ..., verbose = 0,
                          miniter = 5, convergence = TRUE) {
  if (class(x) == "list") {
    x <- array(unlist(x),
      dim = c(
        nrow(x[[1]]),
        ncol(x[[1]]), length(x)
      )
    )
  }
  if (is.null(dim(x))) {
    stop("'x' is not an array")
  }
  if (any(!is.finite(x))) {
    stop("infinite, NA or NaN values in 'x'")
  }
  if (is.null(nu) || nu == 0 || is.infinite(nu)) model <- "normal"

  if (model == "normal") nu <- 0
  # if (model != "normal") {
  #    df = nu
  #
  # }
  df <- nu
  dims <- dim(x)
  ## x is a p * q * n array
  n <- dims[3]
  p <- dims[1]
  q <- dims[2]
  if (verbose > 0) cat("Dims: ", dims, "\n")
  if (!is.null(prior)) {
    if ((length(prior) == 1) && (round(prior) == prior)) {
      prior <- rep(1, prior) / prior
    }

    if (any(prior < 0) || round(sum(prior), 5) != 1) {
      stop("invalid 'prior'")
    }
    prior <- prior[prior > 0L]
    K <- length(prior)
  } else {
    if (missing(K)) stop("No prior and no K")

    prior <- rep(1, K) / K
  }
  if (is.null(init)) {
    init <- init_matrixmixture(x, prior = prior, ...)
  }

  ### extract initialization state
  ### should perhaps handle this by passing to init
  nclass <- length(prior)
  ## if (model != "normal") {
  ## if df is not a vector of length K, take first element and fill out vec
  ## works for normal, too
  if (length(df) != nclass) df <- rep(df[1], nclass)


  centers <- init$centers
  if (!is.null(init$U)) {
    fit_u <- init$U
  } else {
    fit_u <- array(rep(diag(p), nclass), c(p, p, nclass))
    if (model == "t") fit_u <- (df[1] - 2) * stats::var(x[1, 1, ]) * fit_u
  }
  if (!is.null(init$V)) {
    fit_v <- init$V
  } else {
    fit_v <- array(rep(diag(q), nclass), c(q, q, nclass))
  }
  posterior <- matrix(rep(prior, n), byrow = TRUE, nrow = n)
  newposterior <- posterior
  eps <- 1e40
  pi <- prior
  log_lik_vec <- numeric(0)
  if (verbose > 1) {
    cat("\nInit centers: \n\n")
    print(init$centers)
  }
  if (verbose > 2) {
    print("Initial U and V")
    print(fit_u)
    print(fit_v)
  }
  convergeflag <- FALSE

  ss <- array(0, c(p, p, nclass))
  ssx <- array(0, c(p, q, nclass))
  ssxx <- array(0, c(q, q, nclass))
  ssd <- rep(0, nclass)
  new_u <- fit_u
  new_v <- fit_v
  new_df <- df
  newcenters <- centers
  log_lik <- 0
  oldlog_lik <- 0
  olderlog_lik <- 0
  i <- 0
  while (i < iter && (((eps) > tolerance) || (i < miniter))) {
    if (verbose) cat("\nEntering iteration:", i)
    if (verbose > 1) print(pi)
    centers <- newcenters
    newcenters <- array(0, dim = c(p, q, nclass))
    fit_u <- new_u
    fit_v <- new_v
    df <- new_df
    posterior <- newposterior

    ####### E STEP
    ## update expectations of sufficient statistics
    ## update z_ig weights

    for (j in 1:nclass) {
      if (model == "normal") {
        newposterior[, j] <- log(pi[j]) +
          dmatnorm_calc(
            x = x, mean = centers[, , j],
            U = fit_u[, , j], V = fit_v[, , j]
          )
      } else {
        newposterior[, j] <- log(pi[j]) +
          dmat_t_calc(
            x = x,
            df = df[j], mean = centers[, , j],
            U = fit_u[, , j], V = fit_v[, , j]
          )
      }
    }

    newposterior <- ((newposterior - apply(newposterior, 1L,
      min,
      na.rm = TRUE
    )))
    newposterior <- exp(newposterior)
    totalpost <- rowSums(newposterior)
    newposterior <- newposterior / totalpost
    if (verbose > 1) print(newposterior[1:5, ])

    ## update S_ig - conditional weights, only if non-normal

    if (model == "t") {
      dfmult <- df + p + q - 1
      for (j in 1:nclass) {
        s_list <- .sstep(
          x, centers[, , j], fit_u[, , j], fit_v[, , j],
          newposterior[, j]
        )
        ss[, , j] <- s_list$ss
        ssx[, , j] <- s_list$ssx
        ssxx[, , j] <- s_list$ssxx
        ssd[j] <- s_list$ssd
      }
    }

    ### leave blank for now

    ####### CM STEPS
    pi <- colMeans(newposterior)
    if (verbose) cat("\nNew pi: ", pi, "\n")
    ## max for centers, U, V
    ### max for centers

    sumzig <- colSums(newposterior)
    if (verbose > 1) cat("\n Column sums of posterior", sumzig)
    for (j in 1:nclass) {
      newcenters[, , j] <- .means_function(
        x, fit_v[, , j], ss[, , j], ssx[, , j],
        newposterior[, j], row.mean, col.mean, model
      )
    }

    ### max for U, V
    ## if normal
    if (model == "normal") {
      for (j in 1:nclass) {
        ### .NormVarFunc(data,centers,U,V,weights,row.variance,col.variance)
        #### or do EEE, etc formulation
        zigmult <- rep(newposterior[, j], each = q * q)
        swept_data <- sweep(x, c(1, 2), newcenters[, , j])
        inter_v <- txax(swept_data, fit_u[, , j]) * zigmult
        new_v[, , j] <- rowSums(inter_v, dims = 2) / (sumzig[j] * p)
        zigmult <- rep(newposterior[, j], each = p * p)
        inter_u <- xatx(swept_data, new_v[, , j]) * zigmult
        newu <- rowSums(inter_u, dims = 2) / (sumzig[j] * q)
        new_u[, , j] <- newu / (newu[1, 1])
      }
    } else {
      for (j in 1:nclass) {
        new_v[, , j] <- .col_vars(
          x, newcenters[, , j], df[j], newposterior[, j],
          ss[, , j], ssx[, , j], ssxx[, , j], ...
        )$V

        new_u[, , j] <- .row_vars(
          x, newcenters[, , j], df[j], newposterior[, j],
          ss[, , j], ssx[, , j], ssxx[, , j], ...
        )$U
        new_uinv <- (dfmult[j] / (sumzig[j] * (df[j] + p - 1))) * ss[, , j]
        new_u[, , j] <- solve(new_uinv)
      }
    }



    ### Fit NU:
    ### doesn't work yet
    new_df <- df
    ###       if (model == "t" && fixdf == FALSE && iter > 1) {
    ###           ######## THIS DOES NOT WORK.
    ###           for (j in 1:nclass) {
    ###               detss = determinant(ss[,,j], logarithm = TRUE)$modulus[1]
    ###               nu_ll = function(nus) {(CholWishart::mvdigamma((nus + p - 1)/2, p) -
    ###                                     CholWishart::mvdigamma((nus + p + q - 1)/2, p) -
    ###                                    # (ssd[j]/sumzig[j] - (detss - p*log(sumzig[j]*(nus + p - 1))+p*log(nus + p + q - 1))))
    ###                                       # this latest ECME-ish one gives SLIGHTLY different results but is faster
    ###                                       (ssd[j]/sumzig[j] +  determinant(new_u[,,j], logarithm = TRUE)$modulus[1]))
    ###
    ###               }
    ###               if (!isTRUE(sign(nu_ll(2)) * sign(nu_ll(1000)) <= 0)) {
    ###                   warning("Endpoints of derivative of df likelihood do not have opposite sign. Check df specification.")
    ###                   varflag = TRUE
    ###                   ## print(nu_ll(3))
    ###                   ## print(ssd[j])
    ###
    ###               } else {
    ###                   fit0 <- stats::uniroot(nu_ll, c(2, 1000),...)
    ###                   new_df[j] = fit0$root
    ###               }
    ###
    ###           }
    ###
    ###       }

    ####### Eval convergence
    if (verbose > 1) {
      print("New centers:")
      print(newcenters)
      print("New U:")
      print(new_u)
      print("New V:")
      print(new_v)
    }

    olderlog_lik <- oldlog_lik
    oldlog_lik <- log_lik
    log_lik <- 0
    # for(obs in 1:n) {
    for (j in 1:nclass) {
      if (model == "normal" || new_df[j] == 0 || new_df[j] == Inf) {
        log_lik <- log_lik + sum(newposterior[, j] * (log(pi[j]) +
          newposterior[, j] * dmatnorm_calc(
            x = x, mean = newcenters[, , j],
            U = new_u[, , j], V = new_v[, , j]
          )))
      } else {
        log_lik <- log_lik + sum(newposterior[, j] * (log(pi[j]) +
          newposterior[, j] * dmat_t_calc(
            x = x, df = new_df[j], mean = newcenters[, , j],
            U = new_u[, , j], V = new_v[, , j]
          )))
      }
    }
    # }
    if (verbose) cat("\nLog likelihood:", log_lik)
    if (i == 0) {
      oldlog_lik <- log_lik - .3 * abs(log_lik)
      ## initialize to some not-so-bad values
      ## so that doesn't immediately "converge"
      olderlog_lik <- oldlog_lik - .2 * abs(oldlog_lik)
    }

    if (convergence) {
      aitken <- (log_lik - oldlog_lik) / (oldlog_lik - olderlog_lik)
      linf <- oldlog_lik + 1 / (1 - aitken) * (log_lik - oldlog_lik)
      eps <- linf - log_lik
      if (verbose) cat("\nAitken, l_infinity, epsilon:", aitken, linf, eps)
    } else {
      eps <- log_lik - oldlog_lik
    }

    i <- i + 1

    log_lik_vec <- c(log_lik_vec, log_lik)
  }
  if ((i == iter || eps > tolerance)) {
    warning("failed to converge")
  } else {
    convergeflag <- TRUE
  }
  if (verbose) cat("\nDone at iteration ", i - 1, "\n")
  fit_u <- new_u
  fit_v <- new_v
  centers <- newcenters
  posterior <- newposterior
  pi <- colMeans(posterior)
  df <- new_df
  if (verbose > 1) {
    print("Final centers:")
    print(centers)
  }
  if (verbose) cat("\nLog Likelihood Trace: \n", log_lik_vec, "\n")
  cl <- match.call()
  cl[[1L]] <- as.name("matrixmixture")

  structure(
    list(
      prior = prior,
      init = init,
      K = nclass,
      N = n,
      centers = centers,
      U = fit_u,
      V = fit_v,
      posterior = posterior,
      pi = pi,
      nu = df,
      convergence = convergeflag,
      iter = i,
      logLik = log_lik_vec,
      model = model,
      method = method,
      call = cl
    ),
    class = "MixMatrixModel"
  )
}

#' @export
print.MixMatrixModel <- function(x, ...) {
  x[["posterior"]] <- head(x[["posterior"]])
  x[["init"]] <- NULL
  print.default(x, ...)
}

#' @export
#' @importFrom graphics plot
plot.MixMatrixModel <- function(x, ...) {
  plot(
    x = seq_len(length(x$logLik)), y = x$logLik,
    ylab = "Log Likelihood", xlab = "iteration", ...
  )
}


#' @export
logLik.MixMatrixModel <- function(object, ...) {
  dims <- dim(object$centers)
  n <- object$call$N
  p <- dims[1]
  q <- dims[2]
  numgroups <- length(levels(grouping))

  meanpars <- p * q
  if (!is.null(object$call$row.mean) && (object$call$row.mean)) {
    meanpars <- meanpars / q
  }
  if (!is.null(object$call$col.mean) && (object$call$col.mean)) {
    meanpars <- meanpars / p
  }
  upars <- (p + 1) * p / 2
  vpars <- (q + 1) * q / 2 # there's one par that will get subbed off variance
  nupar <- 0 # only fixed for now
  numgroups <- (object$K)
  if (!is.null(object$call$fixdf) && !(object$call$fixdf)) nupar <- numgroups

  ### insert here logic for parsing out different values for this later
  ### as ways of restricting variances and means are added

  df <- numgroups * (vpars + upars + meanpars - 1) + nupar
  log_lik <- object$logLik[length(object$logLik)]

  class(log_lik) <- "logLik"
  attr(log_lik, "df") <- df
  attr(log_lik, "nobs") <- n
  log_lik
}

#' @export
nobs.MixMatrixModel <- function(object, ...) {
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
##' @param prior prior probability. One of \code{prior} and \code{K}
##'      must be  provided. They must be consistent if both provided.
##' @param K number of groups
##' @param centers (optional) either a matrix or an array of
##'      \eqn{p \times p}{p * p}
##'      matrices for use as the \code{centers} argument.
##'      If fewer than \code{K} are provided, the
##'      remainder are chosen by \code{centermethod}.
##' @param U (optional) either a matrix or an array of
##'      \eqn{p \times p}{p * p} matrices for use as the \code{U}
##'      argument. If a matrix is provided, it is duplicated to provide an
##'      array. If an array is provided, it should have \code{K} slices.
##' @param V  (optional) either a matrix or an array of matrices
##'      for use as the \code{V} argument. If a matrix is provided,
##'      it is duplicated to provide an array.
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
##' @param init (optional) a (possibly partially-formed) list
##' with some of the components
##'     \code{centers}, \code{U}, and \code{V}. The function will complete the
##'     list and fill out missing entries.
##' @param ... Additional arguments to pass to \code{kmeans()} if that is
##'     \code{centermethod}.
##' @return a list suitable to use as the \code{init} argument in
##'      \code{matrixmixture}:
##' \describe{
#'       \item{\code{centers}}{the group means,
#'             a \eqn{p \times q \times K}{p * q * K} array.}
#'       \item{\code{U}}{the between-row covariance matrices, a
#'       \eqn{p \times p \times K}{p * p * K}  array}
#'       \item{\code{V}}{the between-column covariance matrix, a
#'       \eqn{q \times q \times K}{q * q * K} array}
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
init_matrixmixture <- function(data, prior = NULL, K = length(prior),
                               centers = NULL, U = NULL, V = NULL,
                               centermethod = "kmeans", varmethod = "identity",
                               model = "normal", init = NULL, ...) {
  dims <- dim(data)
  p <- dims[1]
  q <- dims[2]
  n <- dims[3]
  if (is.null(prior) && is.null(K)) stop("No prior and no K")
  remains <- K
  cenflag <- FALSE
  if (!is.null(centers)) {
    cenflag <- TRUE
    initcenters <- centers
  }
  newcenters <- array(dim = c(p, q, K))
  if (length(prior) == 1) K <- prior
  if (!is.null(K)) prior <- rep(1, K) / K
  if (!is.null(init)) {
    if (!is.null(init$centers)) {
      cenflag <- TRUE
      initcenters <- init$centers
    }
    if (is.null(U)) U <- init$U
    if (is.null(V)) V <- init$V
  }
  if (cenflag) {
    dimcen <- dim(initcenters)
    if (!((dimcen[1] == p) && (dimcen[2] == q))) {
      stop("wrong dimension for provided centers")
    }
    if (length(dimcen) == 2) {
      remains <- K - 1
    } else {
      remains <- K - dimcen[3]
    }
    newcenters[, , (remains + 1):K] <- initcenters
  }
  if (centermethod == "random" && (remains > 0)) {
    select <- sample(n, remains, replace = FALSE)
    newcenters[, , seq_len(remains)] <- data[, , select]
  }
  if ((remains > 0) && (centermethod == "kmeans" ||
    centermethod == "k-means")) {
    res <- kmeans(matrix(data, nrow = n), centers = remains, ...)
    newcenters <- array(res$centers, dim = c(p, q, remains))
    if (cenflag) {
      newcenters <- array(c(newcenters, initcenters),
        dim = c(p, q, K)
      )
    }
  }
  if (!is.null(U)) {
    if (length(dim(U) == 2)) U <- array(rep(U, K), dim = c(p, p, K))
  }
  if (!is.null(V)) {
    if (length(dim(V) == 2)) V <- array(rep(V, K), dim = c(q, q, K))
  }
  if (varmethod == "identity") {
    if (is.null(U)) U <- array(c(rep(diag(p), K)), dim = c(p, p, K))
    if (is.null(V)) V <- array(c(rep(diag(q), K)), dim = c(q, q, K))
  }

  list(
    centers = newcenters,
    U = U,
    V = V
  )
}


# S3 method for predict on class MixMatrixModel
##' @export
predict.MixMatrixModel <- function(object, newdata, prior = object$pi, ...) {
  if (!inherits(object, "MixMatrixModel")) {
    stop("object not of class \"MixMatrixModel\"")
  }

  if (missing(newdata)) {
    newdata <- eval.parent(object$call$x)
  }

  if (is.null(dim(newdata))) {
    stop("'newdata' is not an array")
  }
  if (any(!is.finite(newdata))) {
    stop("infinite, NA or NaN values in 'newdata'")
  }
  x <- (newdata)


  if (length(dim(x)) == 2) x <- array(x, dim = c(dim(x), 1))

  if (ncol(x[, , 1, drop = FALSE]) !=
    ncol(object$centers[, , 1, drop = FALSE])) {
    stop("wrong column dimension of matrices")
  }
  if (nrow(x[, , 1, drop = FALSE]) !=
    nrow(object$centers[, , 1, drop = FALSE])) {
    stop("wrong row dimension of matrices")
  }
  ng <- length(object$prior)
  if (!missing(prior)) {
    if (length(prior) != ng) stop("invalid prior length")
    if (any(prior < 0) || round(sum(prior), 5) != 1) {
      stop("invalid 'prior'")
    }
  }

  dims <- dim(x)
  # x is a p x q x n array
  n <- dims[3]
  p <- dims[1]
  q <- dims[2]
  df <- object$nu


  dist <- matrix(0, nrow = n, ncol = ng)
  posterior <- matrix(0, nrow = n, ncol = ng)




  for (j in seq(ng)) {
    if (object$model == "normal") {
      dist[, j] <- log(prior[j]) + dmatnorm_calc(x,
        mean = matrix(object$centers[, , j], nrow = p, ncol = q),
        U = matrix(object$U[, , j], nrow = p, ncol = p),
        V = matrix(object$V[, , j], nrow = q, ncol = q)
      )
    } else {
      dist[, j] <- log(prior[j]) + dmat_t_calc(x,
        df = df[j],
        mean = matrix(object$centers[, , j], nrow = p, ncol = q),
        U = matrix(object$U[, , j], nrow = p, ncol = p),
        V = matrix(object$V[, , j], nrow = q, ncol = q)
      )
    }
  }
  posterior <- exp((dist - apply(dist, 1L, max, na.rm = TRUE)))
  totalpost <- rowSums(posterior)
  posterior <- posterior / totalpost
  nm <- names(object$prior)
  cl <- max.col(posterior)
  list(class = cl, posterior = posterior)
}
