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
#'     Only positive correlations are allowed for AR(1) and CS and these
#'     restrictions may not be guaranteed to converge.
#'     Note that while maximum likelihood estimators are available (and used)
#'     for the unconstrained variance matrices, \code{optim} is used for any
#'     constraints so it may be considerably slower.
#' @param col.variance  Imposes a variance structure on the columns.
#'     Either 'none', 'AR(1)', 'CS', 'Correlation', or 'Independence'.
#'     Only positive correlations are allowed for
#'     AR(1) and CS.
#' @param df Starting value for the degrees of freedom. If \code{fixed = TRUE},
#'     then this is required and not updated. By default, set to 10.
#' @param fixed Whether \code{df} is estimated or fixed.
#'     By default, \code{TRUE}.
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
#'       \item{\code{convergence}}{a convergence flag,
#'       \code{TRUE} if converged.}
#'       \item{\code{call}}{The (matched) function call.}
#'    }
#'
#' @export
#' @seealso \code{\link{rmatrixnorm}}, \code{\link{rmatrixt}},
#' \code{\link{MLmatrixnorm}}
#'
#' @references
#'     Thompson, G Z.  R Maitra, W Q Meeker, A Bastawros (2019),
#'     "Classification with the matrix-variate-t distribution", arXiv
#'     e-prints arXiv:1907.09565 \url{https://arxiv.org/abs/1907.09565}
#'
#'     Dickey, James M. 1967. “Matricvariate Generalizations of the
#'     Multivariate t Distribution and the Inverted Multivariate t
#'     Distribution.” Ann. Math. Statist. 38 (2): 511–18.
#'     \doi{10.1214/aoms/1177698967}
#'
#'     Liu, Chuanhai, and Donald B. Rubin. 1994. “The ECME Algorithm:
#'     A Simple Extension of EM and ECM with Faster Monotone Convergence.”
#'     Biometrika 81 (4): 633–48.
#'           \doi{10.2307/2337067}
#'
#'    Meng, Xiao-Li, and Donald B. Rubin. 1993. “Maximum Likelihood Estimation
#'    via the ECM Algorithm: A General Framework.” Biometrika 80 (2): 267–78.
#'             \doi{10.1093/biomet/80.2.267}
#'
#'     Rubin, D.B. 1983. “Encyclopedia of Statistical Sciences.” In, 4th ed.,
#'       272–5. John Wiley.
#'
#' @examples
#' set.seed(20180202)
#' # drawing from a distribution with specified mean and covariance
#' A <- rmatrixt(
#'   n = 100, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
#'   L = matrix(c(2, 1, 0, .1), nrow = 2), list = TRUE, df = 5
#' )
#' # fitting maximum likelihood estimates
#' results <- MLmatrixt(A, tol = 1e-5, df = 5)
#' print(results)
MLmatrixt <- function(data, row.mean = FALSE, col.mean = FALSE,
                      row.variance = "none", col.variance = "none",
                      df = 10, fixed = TRUE,
                      tol = .Machine$double.eps^0.5, max.iter = 5000, U, V,
                      ...) {
  if (is.null(df) || df == 0 || is.infinite(df)) {
    return(MLmatrixnorm(
      data, row.mean, col.mean, row.variance,
      col.variance, tol, max.iter, U, V, ...
    ))
  }
  if (class(data) == "list") {
    data <- array(unlist(data),
      dim = c(
        nrow(data[[1]]),
        ncol(data[[1]]), length(data)
      )
    )
  }
  if (!all(
    is.numeric(data), is.numeric(tol),
    is.numeric(max.iter)
  )) {
    stop("Non-numeric input. ")
  }
  if (!(missing(U))) {
    if (!(is.numeric(U))) stop("Non-numeric input.")
  }
  if (!(missing(V))) {
    if (!(is.numeric(V))) stop("Non-numeric input.")
  }

  row_set_var <- FALSE

  rowvarparse <- .varparse(row.variance)
  row_set_var <- rowvarparse$varflag
  row.variance <- rowvarparse$varopt

  col_set_var <- FALSE
  colvarparse <- .varparse(col.variance)
  col_set_var <- colvarparse$varflag
  col.variance <- colvarparse$varopt


  dims <- dim(data)

  if (max(dims[1] / dims[2], dims[2] / dims[1]) > (dims[3] - 1)) {
    stop("Need more observations to estimate parameters.")
  }
  # don't have initial starting point for U and V, start with diag.
  if (missing(U)) {
    est_u <- diag(dims[1])
  } else {
    est_u <- U
  }
  if (missing(V)) {
    est_v <- diag(dims[2])
  } else {
    est_v <- V
  }

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
  swept_data <- sweep(data, c(1, 2), mu)
  iter <- 0
  error_term <- 1e+40

  if (col_set_var) {
    if (est_v[1, 2] > 0) {
      rho_col <- est_v[1, 2]
    } else {
      inter_v <- txax(swept_data, 0.5 * (est_u + t(est_u)))
      est_v <- rowSums(inter_v, dims = 2) / (dims[3] * dims[1])
      if (col.variance == "AR(1)") {
        est_v <- stats::cov2cor(est_v)
        rho_col <- est_v[1, 2]
      }
      if (col.variance == "CS") {
        est_v <- stats::cov2cor(est_v)
        rho_col <- mean(est_v[1, ] / est_v[1, 1])
      }
      if (col.variance == "I") rho_col <- 0
      if (rho_col > .9) rho_col <- .9
      if (rho_col < 0) rho_col <- 0
      est_v <- varmatgenerate(dims[2], rho_col, col.variance)
    }
  }

  if (row_set_var) {
    if (est_u[1, 2] > 0) {
      rho_row <- est_u[1, 2]
    } else {
      inter_u <- xatx(swept_data, 0.5 * (est_v + t(est_v)))
      est_u <- rowSums(inter_u, dims = 2) / (dims[3] * dims[2])
      if (row.variance == "AR(1)") {
        est_u <- stats::cov2cor(est_u)
        rho_row <- est_u[1, 2] / est_u[1, 1]
      }
      if (row.variance == "CS") {
        est_u <- stats::cov2cor(est_u)
        rho_row <- mean(est_u[1, ] / est_u[1, 1])
      }
      if (row.variance == "I") rho_row <- 0
      if (rho_row > .9) rho_row <- .9
      if (rho_row < 0) rho_row <- 0
      est_u <- varmatgenerate(dims[1], rho_row, row.variance)
    }
  }

  varflag <- FALSE

  p <- dims[1]
  q <- dims[2]
  n <- dims[3]


  while (iter < max.iter && error_term > tol && (!varflag)) {
    ### E step
    s_list <- .sstep(data, mu, est_u, est_v, rep(1.0, n))
    ss <- s_list$ss
    ssx <- s_list$ssx
    ssxx <- s_list$ssxx
    ssd <- s_list$ssd


    ### CM STEP
    ### MEANS:
    new_mu <- .means_function(data,
      v = est_v, ss, ssx, rep(1.0, n),
      row.mean, col.mean, "t"
    )

    ### VARS:
    colvarlist <- .col_vars(
      data, new_mu, df, rep(1.0, n), ss, ssx, ssxx,
      col.variance, col_set_var, varflag
    )
    varflag <- colvarlist$varflag
    if (!varflag) new_v <- colvarlist$V

    rowvarlist <- .row_vars(
      data, new_mu, df, rep(1, n), ss, ssx, ssxx,
      row.variance, row_set_var, varflag
    )
    varflag <- rowvarlist$varflag
    if (!varflag) new_u <- rowvarlist$U

    ### IF NU UPDATE
    if (!fixed) {
      new_df <- df
      ## insert E step for NU and M step for NU
      ssdtmp <- ssd

      detss <- determinant(ss, logarithm = TRUE)$modulus[1]
      nu_ll <- function(nu) {
        (CholWishart::mvdigamma((nu + p - 1) / 2, p) -
          CholWishart::mvdigamma((nu + p + q - 1) / 2, p) -
          (ssdtmp / n - (detss - p * log(n * (nu + p - 1)) +
            p * log(nu + p + q - 1))))
        # this latest ECME-ish one gives SLIGHTLY different results but is faster
        # (ssdtmp/n +  determinant(new_u, logarithm = TRUE)$modulus[1]))
      }
      if (!isTRUE(sign(nu_ll(1e-6)) * sign(nu_ll(1000)) <= 0)) {
        warning("Endpoints of derivative of df likelihood do not have
opposite sign.
                 Check df specification.")
        varflag <- TRUE
      } else {
        fit0 <- stats::uniroot(nu_ll, c(1e-6, 1000), ...)
        new_df <- fit0$root
      }
    } else {
      new_df <- df
    }
    ### CHECK CONVERGENCE
    error_term <- sum((new_v - est_v)^2) / (q * q) +
      sum((new_u - est_u)^2) / (p * p) +
      sum((new_mu - mu)^2) / (p * q) + (df - new_df)^2 / (n * p * q)
    ### check, force symmetry
    if (max(abs(new_v - t(new_v)) > tol)) {
      warning("V matrix may not be symmetric")
    }

    if (max(abs(new_u - t(new_u)) > tol)) {
      warning("U matrix may not be symmetric")
    }
    est_v <- .5 * (new_v + t(new_v))
    est_u <- .5 * (new_u + t(new_u))
    mu <- new_mu
    df <- new_df

    iter <- iter + 1
  }
  if (iter >= max.iter || error_term > tol || varflag) {
    warning("Failed to converge")
  }
  log_lik <- sum(dmatrixt(data, mu, U = est_u, V = est_v, df = df, log = TRUE))
  converged <- !(iter >= max.iter || error_term > tol || varflag)

  return(list(
    mean = mu,
    U = est_u / est_u[1, 1],
    V = est_v,
    var = est_u[1, 1],
    nu = df,
    iter = iter,
    tol = error_term,
    logLik = log_lik,
    convergence = converged,
    call = match.call()
  ))
}
