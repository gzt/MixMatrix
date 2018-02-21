

## matnorm classification similar to LDA - create an object like LDA object that can be
## used on vectorized version
matrixlda <-  function(x, grouping, prior, tol = 1.0e-4, ...)  {
    if (class(x) == "list") x <- array(unlist(x),
                                             dim = c(nrow(x[[1]]),
                                                     ncol(x[[1]]), length(x)))
    if(is.null(dim(x))) stop("'x' is not an array")
    if(any(!is.finite(x)))
      stop("infinite, NA or NaN values in 'x'")
    dims = dim(x)
    # x is a p x q x n array
    n <- dims[3]
    p <- dims[1]
    q <- dims[2]
    if(n != length(grouping))
      stop("nrow(x) and length(grouping) are different")
    g <- as.factor(grouping)
    lev <- lev1 <- levels(g)
    counts <- as.vector(table(g))
    if(!missing(prior)) {
      if(any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid 'prior'")
      if(length(prior) != nlevels(g)) stop("'prior' is of incorrect length")
      prior <- prior[counts > 0L]
    }
    if(any(counts == 0L)) {
      empty <- lev[counts == 0L]
      warning(sprintf(ngettext(length(empty),
                               "group %s is empty",
                               "groups %s are empty"),
                      paste(empty, collapse = " ")), domain = NA)
      lev1 <- lev[counts > 0L]
      g <- factor(g, levels = lev1)
      counts <- as.vector(table(g))
    }
    proportions <- counts/n
    ng <- length(proportions)
    names(prior) <- names(counts) <- lev1
    # method <- match.arg(method)
    # if(CV && !(method == "moment" || method == "mle"))
    #   stop(gettext("cannot use leave-one-out CV with method %s",
    #                sQuote(method)), domain = NA)

    group.means = array(0,dim=c(p,q,ng))
    for(i in seq(ng)){
      group.means[,,i] = rowMeans(x[,,g==levels(g)[i]],dims=2 )
    }
    swept.group <- array(0,dims)
    for(i in seq(n)){
      swept.group[,,i] <- x[,,i] - group.means[,,as.numeric(g[i])]
        }
    f1 <- sqrt((apply(swept.group,c(1,2),var)))
    if(any(f1 < tol)) {
      const <- format((1L:(p*q)[f1 < tol]))
      stop(sprintf(ngettext(length(const),
                            "variable %s appears to be constant within groups",
                            "variables %s appear to be constant within groups"),
                   paste(const, collapse = " ")),
           domain = NA)
    }

    varfit <- MLmatrixnorm(swept.group)

    cl <- match.call()
    cl[[1L]] <- as.name("matnormlda")
    structure(list(prior = prior, counts = counts, means = group.means,
                   scaling = varfit$var, U = varfit$U, V = varfit$V, lev = lev, #svd = X.s$d[1L:rank],
                   N = n, call = cl) ,
              class = "matrixlda"
              )
}


predict.matrixlda <- function (object, newdata, prior = object$prior, ...) {
  if (!inherits(object, "matnormlda"))
    stop("object not of class \"matnormlda\"")

  # want this function to ony be internal to this function as it's not "hardened"
  mattrace <- function(x) sum(diag(x))

    if (missing(newdata)) {
      if (!is.null(sub <- object$call$subset))
        newdata <- eval.parent(parse(text = paste(deparse(object$call$x,
                                                          backtick = TRUE), "[", deparse(sub, backtick = TRUE),
                                                  ",]")))
      else newdata <- eval.parent(object$call$x)
      if (!is.null(nas <- object$call$na.action))
        newdata <- eval(call(nas, newdata))
    }

   if(is.null(dim(newdata))) stop("'newdata' is not an array")

  if (class(newdata) == "list") x <- array(unlist(newdata),
                                     dim = c(nrow(newdata[[1]]),
                                             ncol(newdata[[1]]), length(newdata)))
    if(any(!is.finite(newdata)))
      stop("infinite, NA or NaN values in 'newdata'")
    x <- (newdata)

    if (ncol(x[,,1]) != ncol(object$means[,,1]))
      stop("wrong column dimension of matrices")
    if (nrow(x[,,1]) != nrow(object$means[,,1]))
      stop("wrong row dimension of matrices")
    ng <- length(object$prior)
    if (!missing(prior)) {
      if (any(prior < 0) || round(sum(prior), 5) != 1)
        stop("invalid 'prior'")
      if (length(prior) != ng)
        stop("'prior' is of incorrect length")
    }

    dims = dim(x)
    # x is a p x q x n array
    n <- dims[3]
    p <- dims[1]
    q <- dims[2]
  dist = matrix(0,nrow=n, ncol = ng)
  posterior = matrix(0,nrow=n, ncol = ng)
  solveV = solve(object$V * object$scaling)
  solveU = solve(object$U)
  VMUM = numeric(ng)
  for (j in seq(ng))
    VMUM[j] = mattrace( (-.5) * solveV %*% t(object$means[,,j]) %*% solveU %*% object$means[,,j])

  for (i in seq(n)) {
    Xi = x[,,i]
    for (j in seq(ng)) {
      dist[i,j] = mattrace(solveV %*% t(object$means[,,j]) %*% solveU %*% Xi) +  VMUM[j] + log(prior[j])
    }
  }
  posterior = exp(dist)
  totalpost = rowSums(posterior)
  posterior = posterior/totalpost
  nm <- names(object$prior)
  cl <- factor(nm[max.col(posterior)], levels = object$lev)
  list(class = cl, posterior = posterior, x = x)
}


matrixqda <-function(x, grouping, prior, tol = 1.0e-4, ...)  {
  if (class(x) == "list") x <- array(unlist(x),
                                     dim = c(nrow(x[[1]]),
                                             ncol(x[[1]]), length(x)))
  if(is.null(dim(x))) stop("'x' is not an array")
  if(any(!is.finite(x)))
    stop("infinite, NA or NaN values in 'x'")
  dims = dim(x)
  # x is a p x q x n array
  n <- dims[3]
  p <- dims[1]
  q <- dims[2]
  if(n != length(grouping))
    stop("nrow(x) and length(grouping) are different")
  g <- as.factor(grouping)
  lev <- lev1 <- levels(g)
  counts <- as.vector(table(g))
  if(!missing(prior)) {
    if(any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid 'prior'")
    if(length(prior) != nlevels(g)) stop("'prior' is of incorrect length")
    prior <- prior[counts > 0L]
  }
  if(any(counts == 0L)) {
    empty <- lev[counts == 0L]
    warning(sprintf(ngettext(length(empty),
                             "group %s is empty",
                             "groups %s are empty"),
                    paste(empty, collapse = " ")), domain = NA)
    lev1 <- lev[counts > 0L]
    g <- factor(g, levels = lev1)
    counts <- as.vector(table(g))
  }
  proportions <- counts/n
  ng <- length(proportions)
  names(prior) <- names(counts) <- lev1
  # method <- match.arg(method)
  # if(CV && !(method == "moment" || method == "mle"))
  #   stop(gettext("cannot use leave-one-out CV with method %s",
  #                sQuote(method)), domain = NA)

  group.means = array(0,dim=c(p,q,ng))
  groupU = array(0,dim=c(p,p,ng))
  groupV = array(0,dim=c(q,q,ng))
  for(i in seq(ng)){
    # hiding this htere: , ...
    mlfit =  MLmatrixnorm(x[,,g==levels(g)[i]],...)
    if(mlfit$convergence == FALSE) warning("ML fit failed for group ", i)

    group.means[,,i] = mlfit$mean
    groupU[,,i] = mlfit$U
    groupV[,,i] = mlfit$V * mlfit$var
  }
  swept.group <- array(0,dims)
  for(i in seq(n)){
    swept.group[,,i] <- x[,,i] - group.means[,,as.numeric(g[i])]
  }
  f1 <- sqrt((apply(swept.group,c(1,2),var)))
  if(any(f1 < tol)) {
    const <- format((1L:(p*q)[f1 < tol]))
    stop(sprintf(ngettext(length(const),
                          "variable %s appears to be constant within groups",
                          "variables %s appear to be constant within groups"),
                 paste(const, collapse = " ")),
         domain = NA)
  }



  cl <- match.call()
  cl[[1L]] <- as.name("matrixqda")
  structure(list(prior = prior, counts = counts, means = group.means,
                 U = groupU, V = groupV, lev = lev,
                 N = n, call = cl) ,
            class = "matrixqda"
  )
}


predict.matrixqda <- function (object, newdata, prior = object$prior, ...) {
  if (!inherits(object, "matrixqda"))
    stop("object not of class \"matrixqda\"")

  # want this function to ony be internal to this function as it's not "hardened"
  mattrace <- function(x) sum(diag(x))

  if (missing(newdata)) {
    if (!is.null(sub <- object$call$subset))
      newdata <- eval.parent(parse(text = paste(deparse(object$call$x,
                                                        backtick = TRUE), "[", deparse(sub, backtick = TRUE),
                                                ",]")))
    else newdata <- eval.parent(object$call$x)
    if (!is.null(nas <- object$call$na.action))
      newdata <- eval(call(nas, newdata))
  }

  if(is.null(dim(newdata))) stop("'newdata' is not an array")

  if (class(newdata) == "list") x <- array(unlist(newdata),
                                           dim = c(nrow(newdata[[1]]),
                                                   ncol(newdata[[1]]), length(newdata)))
  if(any(!is.finite(newdata)))
    stop("infinite, NA or NaN values in 'newdata'")
  x <- (newdata)

  if (ncol(x[,,1]) != ncol(object$means[,,1]))
    stop("wrong column dimension of matrices")
  if (nrow(x[,,1]) != nrow(object$means[,,1]))
    stop("wrong row dimension of matrices")
  ng <- length(object$prior)
  if (!missing(prior)) {
    if (any(prior < 0) || round(sum(prior), 5) != 1)
      stop("invalid 'prior'")
    if (length(prior) != ng)
      stop("'prior' is of incorrect length")
  }

  dims = dim(x)
  # x is a p x q x n array
  n <- dims[3]
  p <- dims[1]
  q <- dims[2]



  ##### Here is where the work needs to be done.
  dist = matrix(0,nrow=n, ncol = ng)
  posterior = matrix(0,nrow=n, ncol = ng)
  solveU = array(dim = c(p,p,n))
  solveV = array(dim = c(q,q,n))
  for(j in seq(ng)){
  solveV[,,j] = solve(object$V[,,j])
  solveU[,,j] = solve(object$U[,,j])
  }
  VMUM = numeric(ng)
  for (j in seq(ng))
    VMUM[j] = mattrace( (-.5) * solveV[,,j] %*% t(object$means[,,j]) %*% solveU[,,j] %*% object$means[,,j])

  detfactor =  numeric(ng)
  for (j in seq(ng))
    detfactor[j] = .5 * (p * log(det(solveU[,,j])) + n * log(det(solveV[,,j])))

  for (i in seq(n)) {
    Xi = x[,,i]
    for (j in seq(ng)) {
      dist[i,j] = mattrace(-.5 * solveV[,,j] %*% t(Xi) %*% solveU[,,j] %*% Xi) +
        mattrace(solveV[,,j] %*% t(object$means[,,j]) %*% solveU[,,j] %*% Xi) +  VMUM[j] + log(prior[j]) +
        detfactor[j]
    }
  }
  posterior = exp(dist)
  totalpost = rowSums(posterior)
  posterior = posterior/totalpost
  nm <- names(object$prior)
  cl <- factor(nm[max.col(posterior)], levels = object$lev)
  list(class = cl, posterior = posterior, x = x)
}

