#   mlmatrixt.R
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



#' Maximum likelihood estimation for matrix variate t distributions
#'
#' For the matrix variate normal distribution, maximum likelihood estimates
#' exist for \eqn{N > max(p/q,q/p)+1} and are unique for \eqn{N > max(p,q)}.
#' The number necessary for the matrix variate t has not been worked out but
#' this is a lower bound. This implements an ECME algorithm to estimate the
#' mean, covariance, and degrees of freedom parameters. An AR(1), compound
#' symmetry, or independence restriction can be proposed for either or both
#' variance matrices. However, if they are inappropriate for the data, they may
#' fail with a warning.
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
#'     'none', 'AR(1)', 'CS' for 'compound symmetry', 'Correlation' for a
#'     correlation matrix, or 'Independence' for
#'     independent and identical variance across the rows.
#'     Only positive correlations are allowed for AR(1) and CS.
#'     Note that while maximum likelihood estimators are available (and used)
#'     for the unconstrained variance matrices, \code{optim} is used for any
#'     constraints so it may be considerably slower.
#' @param col.variance  Imposes a variance structure on the columns.
#'     Either 'none', 'AR(1)', 'CS', 'Correlation', or 'Independence'.
#'     Only positive correlations are allowed for
#'     AR(1) and CS.
#' @param df Starting value for the degrees of freedom. If \code{fixed = TRUE},
#'     then this is required and not updated. By default, set to 10.
#' @param fixed Whether \code{df} is estimated or fixed. By default, \code{TRUE}.
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
#' @return Returns a list with the following elements:
#' \describe{
#'       \item{\code{mean}}{the mean matrix}
#'       \item{\code{U}}{the between-row covariance matrix}
#'       \item{\code{V}}{the between-column covariance matrix}
#'       \item{\code{var}}{the scalar variance parameter
#'            (the first entry of the covariances are restricted to unity)}
#'       \item{\code{nu}}{the degrees of freedom parameter}
#'       \item{\code{iter}}{the number of iterations}
#'       \item{\code{tol}}{the squared difference between iterations of
#'            the variance matrices at the time of stopping}
#'       \item{\code{logLik}}{log likelihood of result.}
#'       \item{\code{convergence}}{a convergence flag, \code{TRUE} if converged.}
#'       \item{\code{call}}{The (matched) function call.}
#'    }
#'
#' @export
#' @seealso \code{\link{rmatrixnorm}}, \code{\link{rmatrixt}}, \code{\link{MLmatrixnorm}}
#'
#' @references
#'     Dickey, James M. 1967. “Matricvariate Generalizations of the Multivariate t
#'        Distribution and the Inverted Multivariate t
#'        Distribution.” Ann. Math. Statist. 38 (2): 511–18. \doi{10.1214/aoms/1177698967}
#'
#'     Liu, Chuanhai, and Donald B. Rubin. 1994. “The ECME Algorithm: A Simple Extension of
#'           EM and ECM with Faster Monotone Convergence.” Biometrika 81 (4): 633–48.
#'           \doi{10.2307/2337067}
#'
#'    Meng, Xiao-Li, and Donald B. Rubin. 1993. “Maximum Likelihood Estimation via the ECM
#'             Algorithm: A General Framework.” Biometrika 80 (2): 267–78.
#'             \doi{10.1093/biomet/80.2.267}
#' 
#'     Rubin, D.B. 1983. “Encyclopedia of Statistical Sciences.” In, 4th ed., 272–5. John Wiley.

#' 
#' @examples
#' set.seed(20180202)
#' A <- rmatrixt(n=100,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
#'    L=matrix(c(2,1,0,.1),nrow=2),list=TRUE, df = 5)
#' results=MLmatrixt(A, tol = 1e-5, df = 5)
#' print(results)
#'
#'
MLmatrixt <- function(data, row.mean = FALSE, col.mean = FALSE,
                         row.variance = "none", col.variance = "none",
                         df = 10, fixed = TRUE,
                         tol = .Machine$double.eps^0.5, max.iter = 5000, U, V,...) {
  if (is.null(df) || df == 0 || is.infinite(df)) return(MLmatrixnorm(data,row.mean,col.mean,row.variance,col.variance,tol,max.iter,U,V,...))
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
  if (length(row.variance) > 1) stop("Invalid input length for variance: ", row.variance)
  rowvarparse <- .varparse(row.variance)
  row.set.var = rowvarparse$varflag
  row.variance = rowvarparse$varopt
  ## row.set.var = FALSE
  ## if (length(row.variance) > 1) stop("Invalid input length for variance: ", row.variance)
  ## if (grepl("^i", x = row.variance,ignore.case = TRUE)) {
  ##   row.set.var = TRUE
  ##   row.variance = "I"
  ## }
 
  ## if (grepl("^cor", x = row.variance,ignore.case = TRUE)) {
 
  ##   row.variance = "cor"
  ## }
  ## if (grepl("^ar", x = row.variance,ignore.case = TRUE)) {
  ##   row.set.var = TRUE
  ##   row.variance = "AR(1)"
  ## }
  ## if (grepl("^cs", x = row.variance,ignore.case = TRUE)) {
  ##   row.set.var = TRUE
  ##   row.variance = "CS"
  ## }
  col.set.var = FALSE
  if (length(col.variance) > 1) stop("Invalid input length for variance: ", col.variance)

  colvarparse <- .varparse(col.variance)
  col.set.var = colvarparse$varflag
  col.variance = colvarparse$varopt

  
  ## if (grepl("^i", x = col.variance, ignore.case = TRUE)) {
  ##   col.set.var = TRUE
  ##   col.variance = "I"
  ## }
  ## if (grepl("^cor", x = col.variance, ignore.case = TRUE)) {
  
  ##   col.variance = "cor"
  ## }
  ## if (grepl("^ar", x = col.variance, ignore.case = TRUE)) {
  ##   col.set.var = TRUE
  ##   col.variance = "AR(1)"
  ## }
  ## if (grepl("^CS", x = col.variance, ignore.case = TRUE)) {
  ##   col.set.var = TRUE
  ##   col.variance = "CS"
  ## }
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
    mu <- matrix(colMeans(mu), nrow = dims[1], ncol = dims[2], byrow = TRUE)
  }
  # if both are true, this makes it so the mean is constant all over
  swept.data <- sweep(data, c(1, 2), mu)
  iter <- 0
  error.term <- 1e+40

  if (col.set.var) {
    if (V[1,2] > 0) {
      rho.col <- V[1,2]
    } else {

      inter.V <- txax(swept.data, 0.5*(U+t(U)))
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

      inter.U <- xatx(swept.data, 0.5*(V+t(V)))
      U = rowSums(inter.U, dims = 2)/(dims[3]*dims[2])
      if (row.variance == "AR(1)") rho.row <- U[1,2]/U[1,1]
      if (row.variance == "CS") rho.row <- mean(U[1,]/U[1,1])
      if (row.variance == "I") rho.row = 0
      if (rho.row > .9) rho.row <- .9
      if (rho.row < 0) rho.row = 0
      U <- varmatgenerate(dims[1],rho.row,row.variance)
    }
  }

  varflag = FALSE
  logLikvec = numeric(0)
p = dims[1]
q = dims[2]
n = dims[3]
#Smatrix = array(0,c(p,p,n))

  while (iter < max.iter && error.term > tol && (!varflag)) {
    
    dfmult = df + p + q - 1

### E step
      Slist = .SStep(data,mu,U,V,rep(1,n))
      SS = Slist$SS
      SSX = Slist$SSX
      SSXX = Slist$SSXX
      SSD = Slist$SSD

 
### CM STEP
      ### MEANS:
      new.Mu = .MeansFunction(data, U=U,V=V, SS, SSX, rep(1,n), row.mean, col.mean, "t")

      ### VARS:
    
    if (col.variance == "I") {
      new.V = diag(dims[2])
    } else if (col.set.var) {

      nLL <- function(theta) {
        vardetmat <- vardet(dims[2], theta, TRUE, col.variance)
        varinvmat <- varinv(dims[2], theta, TRUE, col.variance)
        #SXOX = rowSums(axbt(SSXtmp,varinvmat,data ), dims = 2)
        SXOX = SSX %*% varinvmat %*% t(rowSums(data, dims = 2)) #only need trace to be equal
        return(-n*p*vardetmat + dfmult  * matrixtrace(SXOX+ SS %*% new.Mu %*% varinvmat %*% t(new.Mu) - SSX %*% varinvmat %*% t(new.Mu) - new.Mu %*% varinvmat %*% t(SSX)))
      }
      if (!isTRUE(sign(nLL(0.01)) * sign(nLL(.999)) <= 0)) {
        warning("Endpoints of derivative of likelihood do not have opposite sign. Check variance specification.")
        rho.col = 0
        varflag = TRUE
      } else {
        fit0 <- stats::uniroot(nLL, c(0.01,.999),...)
        rho.col <- fit0$root
      }
      new.V <- varmatgenerate(dims[2], rho.col,col.variance)
    } else {

      new.V = (dfmult / (n * p)) * (SSXX - t(SSX) %*% new.Mu - t(new.Mu) %*% SSX + t(new.Mu) %*% SS %*% new.Mu)
      if (col.variance == "cor") {
        new.V = stats::cov2cor(new.V)
        if (!all(is.finite(new.V))) {
          varflag = TRUE
          new.V = diag(q)
          }
        } else new.V = new.V/new.V[1,1]
    }
    # Fix V to have unit variance on first component

    if (row.variance == "I") {
      new.U = diag(dims[1]) * n * (df + p - 1)*p/matrixtrace(SS * dfmult)
    } else if (row.set.var) {

      nLL <- function(theta) {
        vardetmat <- vardet(dims[1], theta, TRUE, row.variance)
        Sigma = varmatgenerate(dims[1], theta, row.variance)
        var = n * (df + p - 1)*p/matrixtrace(Sigma %*% SS * dfmult)
        varderivative = varderiv(dims[1],theta, row.variance)
        return(var*dfmult*matrixtrace(varderivative %*% SS) + n*(df + p - 1)*vardetmat)
      }
      if (!isTRUE(sign(nLL(0.01)) * sign(nLL(.999)) <= 0)) {
        warning("Endpoints of derivative of likelihood do not have opposite sign. Check variance specification.")
        rho.row = 0
        varflag = TRUE
      } else {
        fit0 <- stats::uniroot(nLL, c(0.01,.998),...)
        rho.row <- fit0$root
      }

      new.U <- varmatgenerate(dims[1], rho.row,row.variance)
      var = n * (df + p - 1)*p/matrixtrace(new.U %*% SS * dfmult)
      new.U = var*new.U
    } else {

    newUinv = (dfmult/(n * (df + p - 1))) * SS
    new.U = solve(newUinv)
    if (row.variance == "cor") {
      vartmp = exp(mean(log(diag(new.U)))) # should be pos def so no problems
      if (!is.finite(vartmp)) {
        vartmp = 1
        varflag = TRUE
        warning("Variance estimate for correlation matrix not positive definite.")
      }
      new.U = vartmp * stats::cov2cor(new.U)
      # this cute trick preserves the determinant of the matrix
      }
    }

    ### IF NU UPDATE
    if (!fixed) {
    new.df = df
    ## insert E step for NU and M step for NU
    SSDtmp = SSD
    #SSDtmp = detsum(Smatrix)
    detSS = determinant(SS, logarithm = TRUE)$modulus[1]
    nuLL = function(nu) {(CholWishart::mvdigamma((nu + p - 1)/2, p) -
                             CholWishart::mvdigamma((nu + p + q - 1)/2, p) -
                             (SSDtmp/n - (detSS - p*log(n*(nu + p - 1))+p*log(nu + p + q - 1))))
                          # this latest ECME-ish one gives SLIGHTLY different results but is faster
                            #(SSDtmp/n +  determinant(new.U, logarithm = TRUE)$modulus[1]))

    }
    if (!isTRUE(sign(nuLL(1 + 1e-6)) * sign(nuLL(1000)) <= 0)) {
      warning("Endpoints of derivative of df likelihood do not have opposite sign. Check df specification.")
      varflag = TRUE
    }else{
    fit0 <- stats::uniroot(nuLL, c(1 + 1e-6, 1000),...)
    new.df = fit0$root
    }
    #print(new.df)
    } else new.df = df
    ### CHECK CONVERGENCE
    error.term <- sum((new.V - V)^2)/(q*q) + sum((new.U - U)^2)/(p*p) + sum((new.Mu - mu)^2)/(p*q) + (df - new.df)^2/(n*p*q)
    ### check, force symmetry
    if (max(abs(new.V - t(new.V)) > tol)) warning("V matrix may not be symmetric")
        
    if (max(abs(new.U - t(new.U)) > tol)) warning("U matrix may not be symmetric")
    V <- .5*(new.V + t(new.V))
    U <- .5*(new.U + t(new.U))
    mu <- new.Mu
    df <- new.df

    iter <- iter + 1
    
    #logLikvec = c(logLikvec, logLik)
  }
  if (iter >= max.iter || error.term > tol || varflag)
    warning("Failed to converge")
  logLik = sum(dmatrixt(data, mu, U = U, V = V, df = df, log = TRUE))
  converged = !(iter >= max.iter || error.term > tol || varflag)
    
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
