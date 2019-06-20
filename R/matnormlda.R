#' LDA for matrix variate distributions
#'
#' Performs linear discriminant analysis on matrix variate data.
#' This works slightly differently from the LDA function in MASS:
#' it does not sphere the data or otherwise normalize it. It presumes
#' equal variance matrices and probabilities are given as if
#' the data are from a matrix variate normal distribution.
#' The estimated variance matrices are weighted by the prior. However,
#' if there are not enough members of a class to estimate a variance,
#' this may be a problem.
#' The function does not take the formula interface.
#'
#' @family matrixlda
#' @family matrix-variate
#' @param x 3-D array or list of matrix data.
#' @param grouping vector
#' @param prior a vector of prior probabilities of the same length
#'    as the number of classes
#' @param tol by default, \code{1e-4}. Tolerance parameter checks
#'    for 0 variance.
#' @param method whether to use the normal distribution (\code{normal}) or the t-distribution (\code{t}).
#'    By default, normal.
#' @param nu If using the t-distribution, the degrees of freedom parameter. By default, 10.
#' @param ... Arguments passed to or from other methods, such
#'    as additional parameters to pass to \code{MLmatrixnorm} (e.g.,
#'    \code{row.mean})
#'
#'
#' @return Returns a list of class \code{matrixlda} containing
#'    the following components:
#'    \describe{
#'       \item{\code{prior}}{the prior probabilities used.}
#'       \item{\code{counts}}{the counts of group membership}
#'       \item{\code{means}}{the group means.}
#'       \item{\code{scaling}}{the scalar variance parameter}
#'       \item{\code{U}}{the between-row covariance matrix}
#'       \item{\code{V}}{the between-column covariance matrix}
#'       \item{\code{lev}}{levels of the grouping factor}
#'       \item{\code{N}}{The number of observations used.}
#'       \item{\code{method}}{The method used.}
#'       \item{\code{nu}}{The degrees of freedom parameter if the t-distribution was used.}
#'       \item{\code{call}}{The (matched) function call.}
#'    }
#'
#' @export
#'
#' @examples
#' set.seed(20180221)
#' A <- rmatrixnorm(30,mean=matrix(0,nrow=3,ncol=4))
#' B <- rmatrixnorm(30,mean=matrix(1,nrow=3,ncol=4))
#' C <- array(c(A,B), dim=c(3,4,60))
#' groups <- c(rep(1,30),rep(2,30))
#' prior <- c(.5,.5)
#' matrixlda(C, groups, prior)
matrixlda <-  function(x, grouping, prior, tol = 1.0e-4, method = "normal", nu = 10,...)  {
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
  dims = dim(x)
  # x is a p x q x n array
  n <- dims[3]
  p <- dims[1]
  q <- dims[2]
  if (n != length(grouping))
    stop("nrow(x) and length(grouping) are different")
  g <- as.factor(grouping)
  lev <- lev1 <- levels(g)
  counts <- as.vector(table(g))
  if (!missing(prior)) {
    if (any(prior < 0) || round(sum(prior), 5) != 1)
      stop("invalid 'prior'")
    if (length(prior) != nlevels(g))
      stop("'prior' is of incorrect length")
    prior <- prior[counts > 0L]
 }
  if (any(counts == 0L)) {
    empty <- lev[counts == 0L]
    warning(sprintf(
      ngettext(length(empty),
               "group %s is empty",
               "groups %s are empty"),
      paste(empty, collapse = " ")
    ), domain = NA)
    lev1 <- lev[counts > 0L]
    g <- factor(g, levels = lev1)
    counts <- as.vector(table(g))
  }
  proportions <- counts / n
  ng <- length(proportions)
  names(prior) <- names(counts) <- lev1

  group.means = array(0, dim = c(p, q, ng))
  for (i in seq(ng)) {
    group.means[, , i] = suppressWarnings(MLmatrixnorm(x[, , g == levels(g)[i], drop = FALSE], max.iter = 1, ...)$mean)
  }
  swept.group <- array(0, dims)
  for (i in seq(n)) {
    swept.group[, , i] <- x[, , i] - group.means[, , as.numeric(g[i])]
  }
  f1 <- sqrt((apply(swept.group, c(1, 2), stats::var)))
  if (any(f1 < tol)) {
    # this should be caught before this in the MLmatrixnorm stage
    const <- format((1L:(p * q)[f1 < tol]))
    stop(sprintf(
      ngettext(
        length(const),
        "variable %s appears to be constant within groups",
        "variables %s appear to be constant within groups"
      ),
      paste(const, collapse = " ")
    ),
    domain = NA)
  }


  if (method == "t"){
    # for method t with df specified
    Uresult = diag(p)
    Vresult = diag(q)
    varresult = 1
    error = 1e6
    itercount = 0
    while (error > 1e-7 && itercount < 1e4) {
      # this loop is somewhat inelegant
      newUresult = matrix(0,p,p)
      newVresult = matrix(0,q,q)
      newvarresult = 0
      for (i in seq(ng)) {
        varfit <- MLmatrixt(x[, , g == levels(g)[i], drop = FALSE], df = nu,
                            U = Uresult, V = Vresult,...)
        group.means[, , i] <- varfit$mean
        newUresult = newUresult + prior[i] * varfit$U
        newVresult = newVresult + prior[i] * varfit$V
        newvarresult = newvarresult + prior[i] * varfit$var
        if (varfit$convergence == FALSE)
          warning("ML fit failed for group ", i)
      }

      error = sum((newUresult - Uresult)^2) + sum((newVresult - Vresult)^2) + (varresult - newvarresult)^2
      Uresult = newUresult
      Vresult = newVresult
      varresult = newvarresult
      itercount = itercount + 1

    }

  } else {
    Uresult = matrix(0,p,p)
    Vresult = matrix(0,q,q)
    varresult = 0
    #for (i in seq(ng)) {
      varfit <- MLmatrixnorm(swept.group,...)
      Uresult =  varfit$U
      Vresult =  varfit$V
      varresult =  varfit$var
    #}
  }
  cl <- match.call()
  cl[[1L]] <- as.name("matrixlda")
  structure(
    list(
      prior = prior,
      counts = counts,
      means = group.means,
      scaling = varresult,
      U = Uresult,
      V = Vresult,
      lev = lev,
      N = n,
      method = method,
      nu = nu,
      call = cl
    ) ,
    class = "matrixlda"
  )
}


#' Classify Matrix Variate Observations by Linear Discrimination
#'
#' Classify matrix variate observations in conjunction with \code{matrixlda}.
#'
#' This function is a method for the generic function \code{predict()} for
#' class "\code{matrixlda}". It can be invoked by calling predict(x) for
#' an object
#' \code{x} of the appropriate class.
#'
#' @family matrixlda
#' @family matrix-variate
#'
#' @param object object of class "matrixlda"
#' @param newdata array or list of new observations to be classified.
#'     If newdata is missing, an attempt will be made to retrieve the
#'     data used to fit the \code{matrixlda} object.
#' @param prior The prior probabilities of the classes, by default the
#'     proportions in the training set or what was set in the call to
#'     \code{matrixlda}.
#' @param ... arguments based from or to other methods
#'
#' @return
#' Returns a list containing
#'    the following components:
#'    \describe{
#'       \item{\code{class}}{The MAP classification (a factor)}
#'       \item{\code{posterior}}{posterior probabilities for the classes}
#'    }
#'
#' @export
#'
#' @examples
#' set.seed(20180221)
#' A <- rmatrixnorm(30,mean=matrix(0,nrow=3,ncol=4))
#' B <- rmatrixnorm(30,mean=matrix(1,nrow=3,ncol=4))
#' C <- array(c(A,B), dim=c(3,4,60))
#' groups <- c(rep(1,30),rep(2,30))
#' prior <- c(.5,.5)
#' D <- matrixlda(C, groups, prior)
#' predict(D)$posterior[1:10,]
#'
## S3 method for class 'matrixlda'
predict.matrixlda <- function(object, newdata, prior = object$prior, ...) {
    if (!inherits(object, "matrixlda"))
      stop("object not of class \"matrixlda\"")

    # want this function to ony be internal to this function as it's not "hardened"
    mattrace <- function(x)
      sum(diag(x))

    if (missing(newdata)) {
      if (!is.null(sub <- object$call$subset))
        newdata <-
          eval.parent(parse(text = paste(
            deparse(object$call$x,
                    backtick = TRUE),
            "[",
            deparse(sub, backtick = TRUE),
            ",]"
          )))
      else
        newdata <- eval.parent(object$call$x)
      if (!is.null(nas <- object$call$na.action))
        newdata <- eval(call(nas, newdata))
    }

    if (is.null(dim(newdata)))
      stop("'newdata' is not an array")
    if (any(!is.finite(newdata)))
      stop("infinite, NA or NaN values in 'newdata'")
    x <- (newdata)
    if (class(newdata) == "list")
      x <- array(unlist(newdata),
                 dim = c(nrow(newdata[[1]]),
                         ncol(newdata[[1]]), length(newdata)))
    if (length(dim(x)) == 2) x <- array(x, dim= c(dim(x),1))


    if (ncol(x[, , 1, drop = FALSE]) != ncol(object$means[, , 1, drop = FALSE]))
      stop("wrong column dimension of matrices")
    if (nrow(x[, , 1, drop = FALSE]) != nrow(object$means[, , 1, drop = FALSE]))
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
    if (object$method == "t") df = object$nu
    dist = matrix(0, nrow = n, ncol = ng)
    posterior = matrix(0, nrow = n, ncol = ng)
    solveV = matrix(solve(object$V * object$scaling),q,q)
    solveU = matrix(solve(object$U),p,p)
    VMUM = vector("list", ng)
    VMU = vector("list", ng)
    for (j in seq(ng)) {
      VMU[[j]] = solveV %*% crossprod(matrix(object$means[, , j],p,q), solveU )
      VMUM[[j]] =  VMU[[j]] %*% object$means[, , j]
    }

    for (i in seq(n)) {
      Xi = matrix(x[, , i],p,q)
      # if (object$method == "t") UXVX = solveV %*% crossprod(Xi,  solveU) %*% (Xi)
      for (j in seq(ng)) {
        if (object$method == "t") {
          dist[i, j] = -.5 * (df + p + q -1) * log(det(diag(q) + solveV %*% t(Xi - object$means[,,j]) %*% solveU %*% ((Xi - object$means[,,j])))) +
                                                log(prior[j])
        } else dist[i, j] = mattrace(VMU[[j]] %*% Xi) +  mattrace(-.5*VMUM[[j]]) + log(prior[j])
      }
    }

    dist <- ((dist - apply(dist, 1L, max, na.rm = TRUE)))
    posterior = exp(dist)
    totalpost = rowSums(posterior)
    posterior = posterior / totalpost
    nm <- names(object$prior)
    cl <- factor(nm[max.col(posterior)], levels = object$lev)
    list(class = cl, posterior = posterior)
  }


#' Quadratic Discriminant Analysis for Matrix Variate Observations
#'
#' See \code{matrixlda}: quadratic discriminant analysis for matrix
#' variate observations.
#'
#' This uses \code{MLmatrixnorm} to find the means and variances.
#'
#' @family matrixlda
#' @family matrix variate
#' @param x 3-D array or list of matrix data.
#' @param grouping vector
#' @param prior a vector of prior probabilities of the same length
#'    as the number of classes
#' @param tol by default, \code{1e-4}. Tolerance parameter checks
#'    for 0 variance.
#' @param method whether to use the normal distribution (\code{normal}) or the t-distribution (\code{t}).
#'    By default, normal.
#' @param nu If using the t-distribution, the degrees of freedom parameter. By default, 10.
#' @param ... Arguments passed to or from other methods, such
#'    as additional parameters to pass to \code{MLmatrixnorm} (e.g.,
#'    \code{row.mean})
#'
#' @return Returns a list of class \code{matrixqda} containing
#'    the following components:
#'    \describe{
#'       \item{\code{prior}}{the prior probabilities used.}
#'       \item{\code{counts}}{the counts of group membership}
#'       \item{\code{means}}{the group means.}
#'       \item{\code{U}}{the between-row covariance matrices}
#'       \item{\code{V}}{the between-column covariance matrices}
#'       \item{\code{lev}}{levels of the grouping factor}
#'       \item{\code{N}}{The number of observations used.}
#'       \item{\code{method}}{The method used.}
#'       \item{\code{nu}}{The degrees of freedom parameter if the t-distribution was used.}
#'       \item{\code{call}}{The (matched) function call.}
#'    }
#'
#' @export
#'
#' @examples
#' #' set.seed(20180221)
#' A <- rmatrixnorm(30,mean=matrix(0,nrow=3,ncol=4))
#' B <- rmatrixnorm(30,mean=matrix(1,nrow=3,ncol=4))
#' C <- array(c(A,B), dim=c(3,4,60))
#' groups <- c(rep(1,30),rep(2,30))
#' prior <- c(.5,.5)
#' D <- matrixqda(C, groups, prior)
matrixqda <- function(x, grouping, prior, tol = 1.0e-4, method = "normal",  nu = 10, ...)  {
  if (class(x) == "list")
    x <- array(unlist(x),
               dim = c(nrow(x[[1]]),
                       ncol(x[[1]]), length(x)))
  if (is.null(dim(x)))
    stop("'x' is not an array")
  if (any(!is.finite(x)))
    stop("infinite, NA or NaN values in 'x'")
  if (nu == 0 || is.infinite(nu)) method = "normal"

  if (method == "normal") df = NULL
  dims = dim(x)
  # x is a p x q x n array
  n <- dims[3]
  p <- dims[1]
  q <- dims[2]
  if (n != length(grouping))
    stop("nrow(x) and length(grouping) are different")
  g <- as.factor(grouping)
  lev <- lev1 <- levels(g)
  counts <- as.vector(table(g))
  if (!missing(prior)) {
    if (any(prior < 0) ||
        round(sum(prior), 5) != 1)
      stop("invalid 'prior'")
    if (length(prior) != nlevels(g))
      stop("'prior' is of incorrect length")
    prior <- prior[counts > 0L]
  }
  if (any(counts == 0L)) {
    empty <- lev[counts == 0L]
    warning(sprintf(
      ngettext(length(empty),
               "group %s is empty",
               "groups %s are empty"),
      paste(empty, collapse = " ")
    ), domain = NA)
    lev1 <- lev[counts > 0L]
    g <- factor(g, levels = lev1)
    counts <- as.vector(table(g))
  }
  proportions <- counts / n
  ng <- length(proportions)
  names(prior) <- names(counts) <- lev1
  if(method == "t"){
    if(length(nu) != ng)
      df = rep_len(nu, ng) # if you mismatch lengths, you will not have a good time
  }
  group.means = array(0, dim = c(p, q, ng))
  groupU = array(0, dim = c(p, p, ng))
  groupV = array(0, dim = c(q, q, ng))
  for (i in seq(ng)) {
    # hiding this there: , ...
    if (method == "t"){
      mlfit =  MLmatrixt(x[, , g == levels(g)[i], drop = FALSE], nu = df[i], ...)
      df[i] = mlfit$nu
    } else{
      mlfit =  MLmatrixnorm(x[, , g == levels(g)[i], drop = FALSE], ...)
    }
    if (mlfit$convergence == FALSE)
      warning("ML fit failed for group ", i)

    group.means[, , i] = mlfit$mean
    groupU[, , i] = mlfit$U
    groupV[, , i] = mlfit$V * mlfit$var
  }
  swept.group <- array(0, dims)
  for (i in seq(n)) {
    swept.group[, , i] <- x[, , i] - group.means[, , as.numeric(g[i])]
  }
  f1 <- sqrt((apply(swept.group, c(1, 2), stats::var)))
  if (any(f1 < tol)) {
    # this should be caught before this in the MLmatrixnorm stage
    const <- format((1L:(p * q)[f1 < tol]))
    stop(sprintf(
      ngettext(
        length(const),
        "variable %s appears to be constant within groups",
        "variables %s appear to be constant within groups"
      ),
      paste(const, collapse = " ")
    ),
    domain = NA)
  }



  cl <- match.call()
  cl[[1L]] <- as.name("matrixqda")
  structure(
    list(
      prior = prior,
      counts = counts,
      means = group.means,
      U = groupU,
      V = groupV,
      lev = lev,
      N = n,
      method = method,
      nu = df,
      call = cl
    ) ,
    class = "matrixqda"
  )
}


#' Classify Matrix Variate Observations by Quadratic Discrimination
#'
#' Classify matrix variate observations in conjunction with \code{matrixqda}.
#'
#' This function is a method for the generic function \code{predict()} for
#' class "\code{matrixqda}". It can be invoked by calling predict(x) for
#' an object
#' \code{x} of the appropriate class.
#'
#' @family matrixlda
#' 
#'
#' @param object object of class "matrixqda"
#' @param newdata array or list of new observations to be classified.
#'     If newdata is missing, an attempt will be made to retrieve the
#'     data used to fit the \code{matrixqda} object.
#' @param prior The prior probabilities of the classes, by default the
#'     proportions in the training set or what was set in the call to
#'     \code{matrixqda}.
#' @param ... arguments based from or to other methods
#'
#' @return
#' Returns a list containing
#'    the following components:
#'    \describe{
#'       \item{\code{class}}{The MAP classification (a factor)}
#'       \item{\code{posterior}}{posterior probabilities for the classes}
#'    }
#' @export
#'
#' @examples
#' set.seed(20180221)
#' A <- rmatrixnorm(30,mean=matrix(0,nrow=3,ncol=4))
#' B <- rmatrixnorm(30,mean=matrix(1,nrow=3,ncol=4))
#' C <- array(c(A,B), dim=c(3,4,60))
#' groups <- c(rep(1,30),rep(2,30))
#' prior <- c(.5,.5)
#' D <- matrixqda(C, groups, prior)
#' predict(D)$posterior[1:10,]
#'
## S3 method for class "matrixqda"
predict.matrixqda <- function(object, newdata, prior = object$prior, ...) {
    if (!inherits(object, "matrixqda"))
      stop("object not of class \"matrixqda\"")

    # want this function to ony be internal to this function as it's not "hardened"
    mattrace <- function(x)
      sum(diag(x))

    if (missing(newdata)) {
      if (!is.null(sub <- object$call$subset))
        newdata <-
          eval.parent(parse(text = paste(
            deparse(object$call$x,
                    backtick = TRUE),
            "[",
            deparse(sub, backtick = TRUE),
            ",]"
          )))
      else
        newdata <- eval.parent(object$call$x)
      if (!is.null(nas <- object$call$na.action))
        newdata <- eval(call(nas, newdata))
    }

    if (is.null(dim(newdata)))
      stop("'newdata' is not an array")
    if (any(!is.finite(newdata)))
      stop("infinite, NA or NaN values in 'newdata'")
    x <- (newdata)
    if (class(newdata) == "list")
      x <- array(unlist(newdata),
                 dim = c(nrow(newdata[[1]]),
                         ncol(newdata[[1]]), length(newdata)))

    if (length(dim(x)) == 2) x <- array(x, dim= c(dim(x),1))

    if (ncol(x[, , 1, drop = FALSE]) != ncol(object$means[, , 1, drop = FALSE]))
      stop("wrong column dimension of matrices")
    if (nrow(x[, , 1, drop = FALSE]) != nrow(object$means[, , 1, drop = FALSE]))
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
    ##### Here is where the work needs to be done.
    dist = matrix(0, nrow = n, ncol = ng)
    posterior = matrix(0, nrow = n, ncol = ng)
    cholU = vector("list", ng)
    cholV = vector("list", ng)
    solveU = vector("list", ng)
    solveV = vector("list", ng)
    for (j in seq(ng)) {
        cholV[[j]] = chol(object$V[, , j])
        cholU[[j]] = chol(object$U[, , j])
        solveV[[j]] = chol2inv(cholV[[j]])
        solveU[[j]] = chol2inv(cholU[[j]])
    }
    VMUM = vector("list",ng)
    detfactor =  numeric(ng)
    VMU = vector("list",ng)
    for (j in seq(ng)) {
      VMU[[j]] = matrix(solveV[[j]] %*% crossprod(matrix(object$means[, , j],p,q), solveU[[j]]),q,p)
      VMUM[[j]] =  VMU[[j]] %*% matrix(object$means[, , j], p, q)
      logdetU = 2*sum(log(diag(cholU[[j]])))
      logdetV = 2*sum(log(diag(cholV[[j]])))
      detfactor[j] = -.5 * (q * logdetU + p * logdetV)
    }

    for (i in seq(n)) {
      Xi = matrix(x[, , i], p, q)
      for (j in seq(ng)) {
        if (object$method == "t"){

          dist[i, j] = -.5* (df[j] + p + q - 1) * log(det(diag(q) + solveV[[j]] %*% t(Xi - object$means[,,j]) %*% solveU[[j]] %*% (Xi - object$means[,,j]))) + log(prior[j]) +
            detfactor[j]
        } else {
        dist[i, j] = mattrace(-.5 * solveV[[j]] %*% crossprod(Xi, solveU[[j]]) %*% Xi) +
          mattrace(VMU[[j]] %*% Xi) -.5 *  mattrace(VMUM[[j]]) + log(prior[j]) +
          detfactor[j]
        }
      }
    }
    posterior = exp( (dist - apply(dist, 1L, max, na.rm = TRUE)))
    totalpost = rowSums(posterior)
    posterior = posterior / totalpost
    nm <- names(object$prior)
    cl <- factor(nm[max.col(posterior)], levels = object$lev)
    list(class = cl, posterior = posterior)
  }
