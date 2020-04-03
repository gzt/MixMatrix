library(MixMatrix)


context("Testing input dimension integrity")

test_that("Testing bad matrix dimension input", {
  a_mat <- diag(3)
  b_mat <- diag(4)

  expect_error(
    rmatrixnorm(n = 1, mean = a_mat, L = b_mat, R = a_mat, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixnorm(n = 1, mean = a_mat, L = a_mat, R = b_mat, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixnorm(n = 1, mean = a_mat, U = b_mat, R = a_mat, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixnorm(n = 1, mean = a_mat, L = a_mat, V = b_mat, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixnorm(n = 1, mean = a_mat, U = a_mat, V = b_mat, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixnorm(n = 1, mean = a_mat, U = a_mat, R = b_mat, list = FALSE),
    "Non-conforming"
  )

  expect_error(
    rmatrixt(n = 1, df = 1, mean = a_mat, L = b_mat, R = a_mat, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixt(n = 1, df = 1, mean = a_mat, L = a_mat, R = b_mat, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixt(n = 1, df = 1, mean = a_mat, U = b_mat, R = a_mat, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixt(n = 1, df = 1, mean = a_mat, L = a_mat, V = b_mat, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixt(n = 1, df = 1, mean = a_mat, U = a_mat, V = b_mat, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixt(n = 1, df = 1, mean = a_mat, U = a_mat, R = b_mat, list = FALSE),
    "Non-conforming"
  )

  expect_error(
      rmatrixinvt(n = 1, df = 1, mean = a_mat,
                  L = b_mat, R = a_mat, list = FALSE),
    "Non-conforming"
  )

  expect_error(
      rmatrixinvt(n = 1, df = 1, mean = a_mat,
                  L = a_mat, R = b_mat, list = FALSE),
    "Non-conforming"
  )
  expect_error(
      rmatrixinvt(n = 1, df = 1, mean = a_mat,
                  U = b_mat, R = a_mat, list = FALSE),
    "Non-conforming"
  )
  expect_error(
      rmatrixinvt(n = 1, df = 1, mean = a_mat,
                  L = a_mat, V = b_mat, list = FALSE),
    "Non-conforming"
  )
  expect_error(
      rmatrixinvt(n = 1, df = 1, mean = a_mat,
                  U = a_mat, V = b_mat, list = FALSE),
    "Non-conforming"
  )
  expect_error(
      rmatrixinvt(n = 1, df = 1, mean = a_mat,
                  U = a_mat, R = b_mat, list = FALSE),
    "Non-conforming"
  )

  expect_error(
    dmatrixnorm(x = a_mat, mean = a_mat, L = b_mat, R = a_mat),
    "Non-conforming"
  )
  expect_error(
    dmatrixnorm(x = a_mat, mean = a_mat, L = a_mat, R = b_mat),
    "Non-conforming"
  )
  expect_error(
    dmatrixnorm(x = a_mat, mean = a_mat, U = b_mat, R = a_mat),
    "Non-conforming"
  )
  expect_error(
    dmatrixnorm(x = a_mat, mean = a_mat, L = a_mat, V = b_mat),
    "Non-conforming"
  )
  expect_error(
    dmatrixnorm(x = a_mat, mean = a_mat, U = a_mat, V = b_mat),
    "Non-conforming"
  )
  expect_error(
    dmatrixnorm(x = a_mat, mean = a_mat, U = a_mat, R = b_mat),
    "Non-conforming"
  )

  expect_error(
    dmatrixt(x = a_mat, df = 1, mean = a_mat, L = b_mat, R = a_mat),
    "Non-conforming"
  )
  expect_error(
    dmatrixt(x = a_mat, df = 1, mean = a_mat, L = a_mat, R = b_mat),
    "Non-conforming"
  )
  expect_error(
    dmatrixt(x = a_mat, df = 1, mean = a_mat, U = b_mat, R = a_mat),
    "Non-conforming"
  )
  expect_error(
    dmatrixt(x = a_mat, df = 1, mean = a_mat, L = a_mat, V = b_mat),
    "Non-conforming"
  )
  expect_error(
    dmatrixt(x = a_mat, df = 1, mean = a_mat, U = a_mat, V = b_mat),
    "Non-conforming"
  )
  expect_error(
    dmatrixt(x = a_mat, df = 1, mean = a_mat, U = a_mat, R = b_mat),
    "Non-conforming"
  )

  expect_error(
    dmatrixinvt(x = a_mat, df = 1, mean = a_mat, L = b_mat, R = a_mat),
    "Non-conforming"
  )
  expect_error(
    dmatrixinvt(x = a_mat, df = 1, mean = a_mat, L = a_mat, R = b_mat),
    "Non-conforming"
  )
  expect_error(
    dmatrixinvt(x = a_mat, df = 1, mean = a_mat, U = b_mat, R = a_mat),
    "Non-conforming"
  )
  expect_error(
    dmatrixinvt(x = a_mat, df = 1, mean = a_mat, L = a_mat, V = b_mat),
    "Non-conforming"
  )
  expect_error(
    dmatrixinvt(x = a_mat, df = 1, mean = a_mat, U = a_mat, V = b_mat),
    "Non-conforming"
  )
  expect_error(
    dmatrixinvt(x = a_mat, df = 1, mean = a_mat, U = a_mat, R = b_mat),
    "Non-conforming"
  )

  expect_error(MLmatrixnorm(rmatrixnorm(n = 100, mean = a_mat),
                            U = b_mat, V = a_mat))
  expect_error(MLmatrixnorm(rmatrixnorm(n = 100, mean = a_mat),
                            U = a_mat, V = b_mat))
  expect_error(MLmatrixt(rmatrixnorm(n = 100, mean = a_mat),
                         df = 0, U = a_mat, V = b_mat))
  expect_error(MLmatrixt(rmatrixnorm(n = 100, mean = a_mat),
                         U = b_mat, V = a_mat))
  expect_error(MLmatrixt(rmatrixnorm(n = 100, mean = a_mat),
                         U = a_mat, V = b_mat))
})


context("Testing Generators")

test_that("Bad input to generators", {
  n <- 0
  rho <- .5
  expect_error(ARgenerate(n, rho), "greater")
  expect_error(CSgenerate(n, rho), "greater")
  rho <- 1.2
  n <- 2
  expect_error(ARgenerate(n, rho), "less")
  expect_error(CSgenerate(n, rho), "less")
  rho <- -.5
  expect_warning(ARgenerate(n, rho), "greater")
  expect_warning(CSgenerate(n, rho), "greater")
  rho <- .9995
  expect_warning(ARgenerate(n, rho), "correlation")
  expect_warning(CSgenerate(n, rho), "correlation")
  rho <- -2.5
  expect_error(ARgenerate(n, rho), "greater")
  expect_error(CSgenerate(n, rho), "greater")
})


context("Testing LDA")

test_that("Testing bad input to LDA", {
  a_mat <- rmatrixnorm(5, mean = matrix(0, nrow = 2, ncol = 2))
  b_mat <- rmatrixnorm(5, mean = matrix(1, nrow = 2, ncol = 2))
  c_mat <- array(c(a_mat, b_mat), dim = c(2, 2, 10))
  d_mat <- array(0, dim = c(2, 2, 5))
  e_mat <- array(c(a_mat, d_mat), dim = c(2, 2, 10))

  groups <- c(rep(1, 5), rep(2, 5))
  groups_empty <- factor(rep("1", 10), levels = c("1", "2"))
  priors <- c(.5, .5)
  expect_error(
    matrixlda(c(c_mat), grouping = c(rep(1, 40), rep(2, 40)), prior = priors),
    "array"
  )
  expect_error(
    matrixlda(c_mat, grouping = c(rep(1, 120), rep(2, 120)), prior = priors),
    "are different"
  )
  expect_error(
    matrixlda(c_mat, grouping = groups, prior = c(.5, .4)),
    "invalid 'prior'"
  )
  expect_error(
    matrixlda(c_mat, grouping = groups, prior = c(.4, .4, .2)),
    "incorrect length"
  )
  expect_warning(
    matrixlda(c_mat, grouping = groups_empty, prior = priors),
    "empty"
  )
})

context("Testing QDA")



test_that("Testing bad input to QDA", {
  a_mat <- rmatrixnorm(5, mean = matrix(0, nrow = 2, ncol = 2))
  b_mat <- rmatrixnorm(5, mean = matrix(1, nrow = 2, ncol = 2))
  c_mat <- array(c(a_mat, b_mat), dim = c(2, 2, 10))
  d_mat <- array(0, dim = c(2, 2, 5))
  e_mat <- array(c(a_mat, d_mat), dim = c(2, 2, 10))

  groups <- c(rep(1, 5), rep(2, 5))
  groups_empty <- factor(rep("1", 10), levels = c("1", "2"))
  priors <- c(.5, .5)


  expect_error(
    matrixqda(c(c_mat), grouping = c(rep(1, 40), rep(2, 40)), prior = priors),
    "array"
  )
  expect_error(
    matrixqda(c_mat, grouping = c(rep(1, 120), rep(2, 120)), prior = priors),
    "are different"
  )
  expect_error(
    matrixqda(c_mat, grouping = groups, prior = c(.5, .4)),
    "invalid 'prior'"
  )
  expect_error(
    matrixqda(c_mat, grouping = groups, prior = c(.4, .4, .2)),
    "incorrect length"
  )
  expect_warning(
    matrixqda(c_mat, grouping = groups_empty, prior = priors),
    "empty"
  )
  expect_error(
    matrixqda(e_mat, grouping = groups, prior = priors)
  )
})

context("Out of bounds")

test_that("Out of bounds numeric input: ", {
  a_mat <- diag(5)
  a_mat[5, 5] <- 0
  expect_error(rmatrixt(0, 1, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = diag(5)
  ))
  expect_error(rmatrixt(1, -1, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = diag(5)
  ))
  expect_error(rmatrixt(1, 1, matrix(0, nrow = 5, ncol = 5),
    U = -diag(5), V = diag(5)
  ))
  expect_error(rmatrixt(1, 1, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = -diag(5)
  ))
  expect_error(rmatrixt(1, 1, matrix(0, nrow = 5, ncol = 5),
    L = a_mat, V = diag(5)
  ))
  expect_error(rmatrixt(1, 1, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), R = a_mat
  ))
  expect_error(rmatrixt(1, 1, matrix(0, nrow = 5, ncol = 5),
    U = a_mat, V = diag(5)
  ))
  expect_error(rmatrixt(1, 1, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = a_mat
  ))

  expect_error(rmatrixinvt(0, 1, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = diag(5)
  ))
  expect_error(rmatrixinvt(1, -1, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = diag(5)
  ))
  expect_error(rmatrixinvt(1, 1, matrix(0, nrow = 5, ncol = 5),
    U = -diag(5), V = diag(5)
  ))
  expect_error(rmatrixinvt(1, 1, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = -diag(5)
  ))
  expect_error(rmatrixinvt(1, 1, matrix(0, nrow = 5, ncol = 5),
    L = a_mat, V = diag(5)
  ))
  expect_error(rmatrixinvt(1, 1, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), R = a_mat
  ))
  expect_error(rmatrixinvt(1, 1, matrix(0, nrow = 5, ncol = 5),
    U = a_mat, V = diag(5)
  ))
  expect_error(rmatrixinvt(1, 1, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = a_mat
  ))

  expect_error(rmatrixnorm(0, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = diag(5)
  ))
  expect_error(rmatrixnorm(1, matrix(0, nrow = 5, ncol = 5),
    U = -diag(5), V = diag(5)
  ))
  expect_error(rmatrixnorm(1, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = -diag(5)
  ))
  expect_error(rmatrixnorm(1, matrix(0, nrow = 5, ncol = 5),
    L = a_mat, V = diag(5)
  ))
  expect_error(rmatrixnorm(1, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), R = a_mat
  ))
  expect_error(rmatrixnorm(1, 1, matrix(0, nrow = 5, ncol = 5),
    L = a_mat, V = diag(5)
  ))
  expect_error(rmatrixnorm(1, 1, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = a_mat
  ))

  expect_error(dmatrixt(1, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = diag(5)
  ))
  expect_error(dmatrixt(1,
    df = 0, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = diag(5)
  ))

  expect_error(dmatrixinvt(1, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = diag(5)
  ))

  expect_error(dmatrixt(
    df = -1, x = matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = diag(5)
  ))
  expect_error(dmatrixt(
    df = 1, x = matrix(0, nrow = 5, ncol = 5),
    U = -diag(5), V = diag(5)
  ))
  expect_error(dmatrixt(
    df = 1, x = matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = -diag(5)
  ))
  expect_error(dmatrixt(
    df = 1, x = matrix(0, nrow = 5, ncol = 5),
    L = a_mat, V = diag(5)
  ))
  expect_error(dmatrixt(
    df = 1, x = matrix(0, nrow = 5, ncol = 5),
    U = diag(5), R = a_mat
  ))
  expect_error(dmatrixt(
    df = 1, x = matrix(0, nrow = 5, ncol = 5),
    U = a_mat, V = diag(5)
  ))
  expect_error(dmatrixt(
    df = 1, x = matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = a_mat
  ))

  expect_error(dmatrixinvt(
    df = -1, x = matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = diag(5)
  ))
  expect_error(dmatrixinvt(
    df = 1, x = matrix(0, nrow = 5, ncol = 5),
    U = -diag(5), V = diag(5)
  ))
  expect_error(dmatrixinvt(
    df = 1, x = matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = -diag(5)
  ))
  expect_error(dmatrixinvt(
    df = 1, x = matrix(0, nrow = 5, ncol = 5),
    L = a_mat, V = diag(5)
  ))
  expect_error(dmatrixinvt(
    df = 1, x = matrix(0, nrow = 5, ncol = 5),
    U = diag(5), R = a_mat
  ))
  expect_error(dmatrixinvt(
    df = 1, x = matrix(0, nrow = 5, ncol = 5),
    U = a_mat, V = diag(5)
  ))
  expect_error(dmatrixinvt(
    df = 1, x = matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = a_mat
  ))

  expect_error(dmatrixnorm(matrix(0, nrow = 5, ncol = 5),
    U = -diag(5), V = diag(5)
  ))
  expect_error(dmatrixnorm(matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = -diag(5)
  ))
  expect_error(dmatrixnorm(matrix(0, nrow = 5, ncol = 5),
    L = a_mat, V = diag(5)
  ))
  expect_error(dmatrixnorm(matrix(0, nrow = 5, ncol = 5),
    U = diag(5), R = a_mat
  ))
  expect_error(dmatrixnorm(1, matrix(0, nrow = 5, ncol = 5),
    L = a_mat, V = diag(5)
  ))
  expect_error(dmatrixnorm(1, matrix(0, nrow = 5, ncol = 5),
    U = diag(5), V = a_mat
  ))
})
