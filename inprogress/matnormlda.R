

## matnorm classification similar to LDA - create an object like LDA object that can be
## used on vectorized version
matnormlda <-
  function(x, grouping, prior, tol = 1.0e-4, ...)
  {
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
    f1 <- sqrt((apply((x - group.means[, , ]),c(1,2),var)))
    if(any(f1 < tol)) {
      const <- format((1L:(p*q)[f1 < tol]))
      stop(sprintf(ngettext(length(const),
                            "variable %s appears to be constant within groups",
                            "variables %s appear to be constant within groups"),
                   paste(const, collapse = " ")),
           domain = NA)
    }
    scaling <- 1/f1

    ### here is where I left off




    ### need fitting, sphering, etc

    cl <- match.call()
    cl[[1L]] <- as.name("matnormlda")
    structure(list(prior = prior, counts = counts, means = group.means,
                   scaling = scaling, lev = lev, svd = X.s$d[1L:rank],
                   N = n, call = cl),
              class = "lda"
              )
  }
