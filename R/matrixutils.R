#   matrixutils.R
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



#' Generate a unit AR(1) covariance matrix
#'
#' @description generate AR(1) correlation matrices
#'
#' @param n number of columns/rows
#' @param rho correlation parameter
#'
#' @return Toeplitz \eqn{n \times n}{n * n} matrix with 1 on the diagonal
#' and \eqn{rho^k} on the other diagonals, where \eqn{k} is distance from the
#' main diagonal.
#'    Used internally but it is useful for generating your own random matrices.
#' @seealso [stats::toeplitz()]
#' @export
#'
#' @examples
#' ARgenerate(6, .9)
ARgenerate <- function(n, rho) {
  if (n <= 1) {
    stop("n must be greater than 1.")
  }
  if (rho >= 1) {
    stop("rho must be a correlation less than 1.")
  }
  if (rho <= -1) {
    stop("rho must be a correlation greater than -1.")
  }
  if (rho < 0) {
    warning("Rho = ", rho, " and should be greater than 0.")
  }
  if (rho > 0.999) {
    warning("Rho = ", rho, " high correlation may cause numerical problems.")
  }
  stats::toeplitz(c(1, rho^(1:(n - 1))))
}


#' Generate a compound symmetric correlation matrix
#'
#' @param n number of dimensions
#' @param rho off-diagonal element - a correlation between -1 and 1.
#' Will warn if less than 0.
#'
#' @return returns an \eqn{n \times n}{n * n} matrix with 1 on the diagonal and
#' `rho` on the off-diagonal.
#' @export
#'
#' @examples
#' # generates a covariance matrix with 1 on the main diagonal
#' # and 0.5 on the off-diagonal elements.
#' CSgenerate(3, .5)
CSgenerate <- function(n, rho) {
  if (n <= 1) {
    stop("n must be greater than 1.")
  }
  if (rho >= 1) {
    stop("rho must be a correlation less than 1.")
  }
  if (rho <= -1) {
    stop("rho must be a correlation greater than -1.")
  }
  if (rho < 0) {
    warning("Rho = ", rho, " and should be greater than 0.")
  }
  if (rho > 0.999) {
    warning("Rho = ", rho, " high correlation may cause numerical problems.")
  }
  mat <- matrix(rho, nrow = n, ncol = n)
  diag(mat) <- 1
  mat
}


#' Quick symmetry check
#'
#' Quick check whether matrix input is symmetric -
#' checks sum of absolute differences of transposes
#' as well as dimensions. Not robust, so only an
#' internal function to be used with known safe input.
#'
#' @param mat Numeric real matrix. Does not check if real.
#' @param tol tolerance - note that if you have a big matrix
#'    it may need to be specified as it's a sum of entries.
#' @noRd
#' @return logical TRUE if symmetric FALSE otherwise.
#'     Not as robust as `isSymmetric()`.
#' @keywords internal
symm_check <- function(mat, tol = (.Machine$double.eps)^.5) {
  if (!is.matrix(mat)) {
    return(FALSE)
  }
  if (!is.numeric(mat)) {
    return(FALSE)
  }
  # commented those out because it is always checked before running symm_check.
  dims <- dim(mat)
  if (dims[1] != dims[2]) {
    return(FALSE)
  }

  testsymmetric(mat, tol)
}


#' Select a variance structure to generate.
#'
#' @param n number of dimensions
#' @param rho parameter for selected variance structure.
#' @param variance variance structure - AR(1) or CS.
#' @noRd
#' @return Specified matrix structure
#'
#' @keywords internal
varmatgenerate <- function(n, rho, variance) {
  if (variance == "I" || variance == "independence" ||
    variance == "Independence") {
    variance <- "I"
  }
  if (variance == "AR(1)") {
    return(ARgenerate(n, rho))
  }
  if (variance == "CS") {
    return(CSgenerate(n, rho))
  }
  if (variance == "I") {
    return(diag(n))
  } else {
    stop("Bad covariance structure input.", variance)
  }
}

#' Determinant selector for chosen covariance matrix.
#'
#' @param n dimensions
#' @param rho off-diagonal parameter
#' @param deriv logical whether to return the determinant or the derivative of
#'     the log of the determinant
#' @param variance  variance structure - AR(1) or CS.
#'
#' @return Determinant or derivative of log-inverse for the
#' specified matrix structure.
#' @keywords internal
vardet <- function(n, rho, deriv, variance) {
  if (variance == "AR(1)") {
    return(ar_det(n, rho, deriv))
  }

  if (variance == "CS") {
    return(cs_det(n, rho, deriv))
  }

  if (variance == "I") {
    return(1)
  } else {
    stop("Bad covariance structure input.", variance)
  }
}

#' Inverse selector for chosen covariance matrix.
#'
#' @param n dimensions
#' @param rho off-diagonal parameter
#' @param deriv logical whether to return the inverse or the derivative of
#'     the inverse
#' @param variance  variance structure - AR(1) or CS.
#' @noRd
#' @return The inverse or derivative of the inverse of the
#' selected matrix structure.
#' @keywords internal
varinv <- function(n, rho, deriv, variance) {
  if (variance == "AR(1)") {
    return(inv_ar(n, rho, deriv))
  }

  if (variance == "CS") {
    return(inv_cs(n, rho, deriv))
  }

  if (variance == "I") {
    return(diag(n))
  } else {
    stop("Bad covariance structure input. ", variance)
  }
}


#' Determinant for an AR(1) covariance matrix.
#'
#' @param n dimensions
#' @param rho off-diagonal parameter
#' @param deriv logical whether to return the determinant or the derivative of
#'     the log of the determinant
#' @noRd
#' @return determinant of an AR(1) covariance matrix.
#'     If `deriv` is specified,
#'     will return the derivative of \eqn{\log |\Sigma^{-1}|}.
#' @keywords internal
ar_det <- function(n, rho, deriv = FALSE) {
  if (!deriv) {
    return((1 - rho^2)^(n - 1))
  } else {
    (1 - n) * (-2 * rho) / (1 - rho^2)
  }
}


#' Determinant for an CS covariance matrix.
#'
#' @param n dimensions
#' @param rho off-diagonal parameter
#' @param deriv logical whether to return the determinant or the
#' derivative of the log determinant of the inverse
#' @noRd
#' @return determinant of an CS covariance matrix
#' @keywords internal
cs_det <- function(n, rho, deriv = FALSE) {
  if (!deriv) {
    return((1 + rho * (n - 1)) * (1 - rho)^(n - 1))
  } else {
    (n - 1) * n * rho / ((1 - rho) * ((n - 1) * rho + 1))
  }
}


#' Returns the inverse of an AR(1) covariance matrix or its derivative
#'
#' @param n dimensions of matrix
#' @param rho correlation parameter
#' @param deriv logical. if TRUE will output the derivative of the
#' inverse matrix.
#' @noRd
#' @return Matrix of the inverse of the AR(1) covariance matrix (or its inverse)
#' @keywords internal
inv_ar <- function(n, rho, deriv = FALSE) {
  if (!(n > 1)) stop("n needs to be greater than 1")
  if (!(rho < 1 && rho > -1)) stop("rho needs to be < 1")
  if (deriv) {
    #-(rho^2+1)/(1-rho^2)^2  , 4*rho/((1-rho^2)^2), 2*rho/((1-rho^2)^2)
    mat <- stats::toeplitz(c(
      4 * rho,
      -(rho^2 + 1),
      rep(0, n - 2)
    ))
    mat[1, 1] <- 2 * rho
    mat[n, n] <- 2 * rho
    return((1 / (1 - rho^2)^2) * mat)
  }
  mat <- stats::toeplitz(c(1 + rho^2, -rho, rep(0, n - 2)))
  mat[1, 1] <- mat[n, n] <- 1
  ((1 / (1 - rho^2)) * mat)
}

#' Inverse of Compound Symmetric Matrix.
#'
#' @param n dimension of the matrix
#' @param rho compound symmetric factor
#' @param deriv logical, whether to return the derivative of the inverse.
#' @noRd
#' @return inverse matrix or derivative of the inverse matrix.
#' @keywords internal
#'
inv_cs <- function(n, rho, deriv = FALSE) {
  if (!(n > 1)) stop("n needs to be greater than 1")
  if (!(rho < 1 && rho > -1)) stop("rho needs to be < 1")
  alpha <- sqrt(1 - rho)

  if (deriv) {
    mat <- diag(1 / alpha^4, n) -
      matrix(((n - 1) * rho^2 + 1) / ((1 - rho)^2 * ((n - 1) * rho + 1)^2),
        nrow = n, ncol = n
      )
    return(mat)
  }

  diag(1 / alpha^2, n) - matrix(rho / (alpha^2 * (alpha^2 + n * rho)),
    nrow = n, ncol = n
  )
}



#' Derivative selector for chosen covariance matrix.
#'
#' @param n dimensions
#' @param rho off-diagonal parameter
#' @param variance  variance structure - AR(1) or CS.
#' @noRd
#' @return The derivative of the selected matrix structure.
#' @keywords internal
varderiv <- function(n, rho, variance) {
  if (variance == "AR(1)") {
    return(deriv_ar(n, rho))
  }

  if (variance == "CS") {
    return(deriv_cs(n, rho))
  }

  if (variance == "I") {
    return(matrix(0, n, n))
  } else {
    stop("Bad covariance structure input. ", variance)
  }
}

#' Derivative of AR(1) covariance matrix
#'
#' @param n dimensions
#' @param rho off-diagonal parameter
#' @return matrix derivative of an AR(1) matrix
#' @noRd
#' @keywords internal
deriv_ar <- function(n, rho) {
  if (!(n > 1)) stop("n needs to be greater than 1")
  if (!(rho < 1 && rho > -1)) stop("rho needs to be < 1")
  x <- 0:(n - 1)
  y <- rho^(-1:(n - 2))
  stats::toeplitz(x * y)
}


#' Derivative of compound symmetric covariance matrix
#'
#' @param n dimensions
#' @param rho off-diagonal parameter
#' @return matrix derivative of a CS matrix with unit diagonal
#' @noRd
#' @keywords internal
deriv_cs <- function(n, rho) {
  if (!(n > 1)) stop("n needs to be greater than 1")
  if (!(rho < 1 && rho > -1)) stop("rho needs to be < 1")

  stats::toeplitz(c(0, rep(1, (n - 1))))
}



#' Trace of a matrix
#'
#' @param mat square numeric matrix
#' @noRd
#' @return trace of a matrix
#' @keywords internal
matrixtrace <- function(mat) {
  b <- dim(mat)
  if (b[1] != b[2]) warning("non-conformable dimensions")
  result <- sum(diag(mat))
  return(result)
}

#' Matrix of ones
#'
#' @param nrow Number of rows
#' @param ncol Number of columns
#' @keywords internal
#' @noRd
#' @return Matrix of ones
ones <- function(nrow, ncol) {
  matrix(1, nrow = nrow, ncol = ncol)
}
