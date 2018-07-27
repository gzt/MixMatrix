
#' MLmatrixt:
#'
#' @description Maximum likelihood estimation for matrix normal distributions
#'
#' Maximum likelihood estimates exist for \eqn{N > max(p/q,q/p)+1} and are
#' unique for \eqn{N > max(p,q)}. This finds the estimate for the mean and then alternates
#' between estimates for the \eqn{U} and \eqn{V} matrices until convergence.
#' An AR(1), compound symmetry, or independence restriction can be proposed for either or both
#' variance matrices. However, if they are inappropriate for the data, they may fail with
#' a warning.
#'
#' @param data Either a list of matrices or a 3-D array with matrices in
#'    dimensions 1 and 2, indexed by dimension 3.
#' @param row.mean By default, \code{FALSE}. If \code{TRUE}, will fit a
#'    common mean within each row. If both this and \code{col.mean} are
#'    \code{TRUE}, there will be a common mean for the entire matrix.
#' @param col.mean By default, \code{FALSE}. If \code{TRUE}, will fit a
#'    common mean within each row. If both this and \code{row.mean} are
#'    \code{TRUE}, there will be a common mean for the entire matrix.
#' @param row.variance Imposes a variance structure on the rows. Either
#'     'none', 'AR(1)', 'CS' for 'compound symmetry', or 'Independence' for
#'     independent and identical variance across the rows.
#'     Only positive correlations are allowed for AR(1) and CS.
#'     Note that while maximum likelihood estimators are available (and used) for
#'     the unconstrained variance matrices, \code{optim} is used for any
#'     constraints so it may be considerably slower.
#' @param col.variance  Imposes a variance structure on the columns.
#'     Either 'none', 'AR(1)', 'CS', or 'Independence'. Only positive correlations are allowed for
#'     AR(1) and CS.
#' @param df Starting value for the degrees of freedom. If fixed = TRUE, then this is required and not updated. By default, set to 10.
#' @param fixed Whether nu is estimated or fixed. By default, TRUE.
#' @param tol Convergence criterion. Measured against square deviation
#'    between iterations of the two variance-covariance matrices.
#' @param max.iter Maximum possible iterations of the algorithm.
#' @param U (optional) Can provide a starting point for the U matrix.
#'    By default, an identity matrix.
#' @param V (optional) Can provide a starting point for the V matrix.
#'    By default, an identity matrix.
#' @param ... (optional) additional arguments can be passed to \code{optim}
#'    if using restrictions on the variance.
#'
#' @return Returns a list with a mean matrix, a \eqn{U} matrix, a \eqn{V}
#'    matrix, the variance parameter (the first entry of the variance matrices
#'    are constrained to be 1 for uniqueness), the number of iterations, the squared difference
#'    between iterations of the variance matrices at the time of stopping, the log likelihood,
#'    and a convergence code.
#' @export
#' @seealso \code{\link{rmatrixnorm}}
#'
#' @examples
#' set.seed(20180202)
#' A <- rmatrixnorm(n=100,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'    L=matrix(c(2,1,0,.1),nrow=2),list=TRUE)
#' results=MLmatrixnorm(A, tol = 1e-5)
#' print(results)
#'
#'
MLmatrixt <- function(data, row.mean = FALSE, col.mean = FALSE,
                         row.variance = "none", col.variance = "none",
                         df = 10, fixed = TRUE,
                         tol = .Machine$double.eps^0.5, max.iter = 5000, U, V,...) {
  if (df == 0 || is.infinite(df)) return(MLmatrixnorm(data,row.mean,col.mean,row.variance,col.variance,tol,max.iter,U,V,...))
    if (class(data) == "list") data <- array(unlist(data),
                                           dim = c(nrow(data[[1]]),
                                                   ncol(data[[1]]), length(data)))
  if (!all(is.numeric(data),is.numeric(tol),
           is.numeric(max.iter))) stop("Non-numeric input. ")
  if (!(missing(U))) {
    if (!(is.numeric(U))) stop("Non-numeric input.")
  }
  if (!(missing(V))) {
    if (!(is.numeric(V))) stop("Non-numeric input.")
  }
  row.set.var = FALSE
  if (length(row.variance) > 1) stop("Invalid input length for variance: ", row.variance)
  if (grepl("^i", x = row.variance,ignore.case = T)) {
    row.set.var = TRUE
    row.variance = "I"
  }
  if (row.variance == "AR(1)" || row.variance == "CS") row.set.var = TRUE

  col.set.var = FALSE
  if (length(col.variance) > 1) stop("Invalid input length for variance: ", col.variance)
  if (grepl("^i", x = col.variance, ignore.case = T)) {
    col.set.var = TRUE
    col.variance = "I"
  }
  if (col.variance == "AR(1)" || col.variance == "CS" ) col.set.var = TRUE
  # if data is array, presumes indexed over third column (same as output
  # of rmatrixnorm) if list, presumes is a list of the matrices
  dims <- dim(data)

  if (max(dims[1]/dims[2], dims[2]/dims[1]) > (dims[3] - 1))
    warning("Need more observations to estimate parameters.")
  # don't have initial starting point for U and V, start with diag.
  if (missing(U))
    U <- diag(dims[1])
  if (missing(V))
    V <- diag(dims[2])
    # mu <- apply(data, c(1, 2), mean)
  mu <- rowMeans(data, dims = 2)
  if (row.mean) {
    # make it so that the mean is constant within a row
    mu <- matrix(rowMeans(mu), nrow = dims[1], ncol = dims[2])
  }
  if (col.mean) {
    # make it so that the mean is constant within a column
    mu <- matrix(colMeans(mu), nrow = dims[1], ncol = dims[2], byrow = T)
  }
  # if both are true, this makes it so the mean is constant all over
  swept.data <- sweep(data, c(1, 2), mu)
  iter <- 0
  error.term <- 1e+40
  if (col.set.var) {
    if (V[1,2] > 0) {
      rho.col <- V[1,2]
    } else {

      inter.V <- txax(swept.data, U)
      V <- rowSums(inter.V, dims = 2)/(dims[3] * dims[1])
      if (col.variance == "AR(1)") rho.col <- V[1,2]/V[1,1]
      if (col.variance == "CS") rho.col <- mean(V[1,]/V[1,1])
      if (col.variance == "I") rho.col = 0
      if (rho.col > .9) rho.col <- .9
      if (rho.col < 0) rho.col <- 0
      V <- varmatgenerate(dims[2],rho.col,col.variance)
    }
  }

  if (row.set.var) {
    if (U[1,2] > 0) {
      rho.row <- U[1,2]
    } else {

      inter.U <- xatx(swept.data, V)
      U = rowSums(inter.U, dims = 2)/(dims[3]*dims[2])
      if (row.variance == "AR(1)") rho.row <- U[1,2]/U[1,1]
      if (row.variance == "CS") rho.row <- mean(U[1,]/U[1,1])
      if (row.variance == "I") rho.row = 0
      if (rho.row > .9) rho.row <- .9
      if (rho.row < 0) rho.row = 0
      U <- varmatgenerate(dims[1],rho.row,row.variance)
    }
  }

matrixtrace = function(A){
  b = dim(A)
  result = 0
  if (b[1] != b[2]) warning("non-conformable dimensions")
    else {
      for (i in 1:b[1]) {
        result = result + A[i,i]
      }
    }
  return(result)
}

  varflag = FALSE
p = dims[1]
q = dims[2]
n = dims[3]
Smatrix = array(0,c(p,p,n))

  while (iter < max.iter && error.term > tol && (!varflag)) {
    swept.data <- sweep(data, c(1, 2), mu)
    dfmult = df + p + q - 1

    ### E step

    Stmp = xatx(swept.data,V)
    for (i in 1:n) Stmp[,,i] = Stmp[,,i] + U
    Smatrix = cubeinv(Stmp)

    SS = rowSums(Smatrix,FALSE, 2)

    SSXtmp = cubemult(Smatrix, data)
    SSX = rowSums(SSXtmp, FALSE, 2)

    SSXXtmp = cubemult(data,SSXtmp)
    SSXX = rowSums(SSXXtmp,FALSE, 2)
    #print(SS)
    #print(SSX)
    #print(SSXX)
    #print(iter)

    ### CM STEP

      if (row.mean && col.mean) {
        # make it so that the mean is constant within a row
        scalarmu = matrixtrace(SSX %*% solve(V) %*% matrix(1, nrow = q, ncol = p)) / matrixtrace(SS %*% matrix(1, nrow = p, ncol = q) %*% solve(V) %*% matrix(1, nrow = q, ncol = p))
        new.Mu <-   scalarmu * matrix(1, nrow = p, ncol = q)
        } else if (col.mean) {
        # make it so that the mean is constant within a column
        # ie mu = p x 1, times ones 1 x q
        new.Mu <- matrix(1, nrow = p, ncol = p) %*% SSX / sum(SS)
        } else if (row.mean) {
          # make it so that the mean is constant within a row
          # ie  ones p x 1 times mu = 1 x q
        new.Mu = solve( SS) %*% SSX %*% (solve(V) %*% matrix(1, nrow = q, ncol = q)) / sum(solve(V))
          } else {
    new.Mu =  solve( SS) %*% SSX
      }
    #new.V = (dfmult / (n * p)) * (SSXX - t(SSX) %*% solve(SS) %*% (SSX))

    if (col.variance == "I") {
      new.V = diag(dims[2])
    } else if (col.set.var) {

      nLL <- function(theta) {
        vardetmat <- vardet(dims[2], theta, TRUE, col.variance)
        varinvmat <- varinv(dims[2], theta, TRUE, col.variance)
        SXOX = rowSums(axbt(SSXtmp,varinvmat,data ),dims=2)

        return(-n*p*vardetmat + dfmult  * matrixtrace(SXOX+ SS %*% new.Mu %*% varinvmat %*% t(new.Mu) - SSX %*% varinvmat %*% t(new.Mu) - new.Mu %*% varinvmat %*% t(SSX)))
      }
      if (!isTRUE(sign(nLL(0)) * sign(nLL(.999)) <= 0)) {
        warning("Endpoints of derivative of likelihood do not have opposite sign. Check variance specification.")
        rho.col = 0
        varflag = TRUE
      } else {
        fit0 <- stats::uniroot(nLL, c(0,.999),...)
        rho.col <- fit0$root
      }
      new.V <- varmatgenerate(dims[2], rho.col,col.variance)
    } else {

      new.V = (dfmult / (n * p)) * (SSXX - t(SSX) %*% new.Mu - t(new.Mu) %*% SSX + t(new.Mu) %*% SS %*% new.Mu)
      new.V = new.V/new.V[1,1]
    }



    # Fix V to have unit variance on first component

    newUinv = (dfmult/(n * (df + p - 1))) * SS
    new.U = solve(newUinv)




    ### IF NU UPDATE
    if (!fixed) {
    new.df = df
    ## insert E step for NU and M step for NU

    SSDtmp = detsum(Smatrix)
    detSS = determinant(SS, logarithm = TRUE)$modulus[1]
    nuLL = function(nu) {(CholWishart::mvdigamma((nu + p - 1)/2, p) -
                             CholWishart::mvdigamma((nu + p + q - 1)/2, p) -
                            (SSDtmp/n - (detSS - p*log(n*(nu + p - 1)/(nu + p + q - 1)))))
                          # this latest ECME-ish one gives SLIGHTLY different results but is faster
                          #  (SSDtmp + n * determinant(new.U)$modulus[1] ))
    }
    if (!isTRUE(sign(nuLL(p - 1)) * sign(nuLL(1000)) <= 0)) {
      warning("Endpoints of derivative of df likelihood do not have opposite sign. Check df specification.")
      varflag = TRUE
    }else{
    fit0 <- stats::uniroot(nuLL, c(p - 1, 1000),...)
    new.df = fit0$root
    }
    #print(new.df)
    } else new.df = df
    ### CHECK CONVERGENCE
    error.term <- sum((new.V - V)^2)/(q*q) + sum((new.U - U)^2)/(p*p) + sum((new.Mu - mu)^2/(p*q) + 1/(n*p*q) * (df - new.df)^2)
    V <- new.V
    U <- new.U
    mu <- new.Mu
    df <- new.df


    iter <- iter + 1
  }
  if (iter >= max.iter || error.term > tol || varflag)
    warning("Failed to converge")

  converged = !(iter >= max.iter || error.term > tol || varflag)
  logLik = 0
  #for (i in seq(dims[3])) {
  #  logLik = logLik + dmatrixnorm(data[,,i], mu, U = U, V = V, log = TRUE)
  #}
  logLik = sum(dmatrixt(data, mu, U = U, V = V, df = df, log = TRUE))
  return(list(mean = mu,
              U = U/U[1,1],
              V = V,
              var = U[1,1],
              nu = df,
              iter = iter,
              tol = error.term,
              logLik = logLik,
              convergence = converged,
              call = match.call()))
}
