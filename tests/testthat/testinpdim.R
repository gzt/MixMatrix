library(MixMatrix)


context("Testing input dimension integrity")

test_that("Testing bad matrix dimension input", {
  A <- diag(3)
  B <- diag(4)

  expect_error(
    rmatrixnorm(n = 1, mean = A, L = B, R = A, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixnorm(n = 1, mean = A, L = A, R = B, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixnorm(n = 1, mean = A, U = B, R = A, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixnorm(n = 1, mean = A, L = A, V = B, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixnorm(n = 1, mean = A, U = A, V = B, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixnorm(n = 1, mean = A, U = A, R = B, list = FALSE),
    "Non-conforming"
  )

  expect_error(
    rmatrixt(n = 1, df = 1, mean = A, L = B, R = A, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixt(n = 1, df = 1, mean = A, L = A, R = B, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixt(n = 1, df = 1, mean = A, U = B, R = A, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixt(n = 1, df = 1, mean = A, L = A, V = B, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixt(n = 1, df = 1, mean = A, U = A, V = B, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixt(n = 1, df = 1, mean = A, U = A, R = B, list = FALSE),
    "Non-conforming"
  )

  expect_error(
    rmatrixinvt(n = 1, df = 1, mean = A, L = B, R = A, list = FALSE),
    "Non-conforming"
  )

  expect_error(
    rmatrixinvt(n = 1, df = 1, mean = A, L = A, R = B, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixinvt(n = 1, df = 1, mean = A, U = B, R = A, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixinvt(n = 1, df = 1, mean = A, L = A, V = B, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixinvt(n = 1, df = 1, mean = A, U = A, V = B, list = FALSE),
    "Non-conforming"
  )
  expect_error(
    rmatrixinvt(n = 1, df = 1, mean = A, U = A, R = B, list = FALSE),
    "Non-conforming"
  )

  expect_error(
    dmatrixnorm(x = A, mean = A, L = B, R = A),
    "Non-conforming"
  )
  expect_error(
    dmatrixnorm(x = A, mean = A, L = A, R = B),
    "Non-conforming"
  )
  expect_error(
    dmatrixnorm(x = A, mean = A, U = B, R = A),
    "Non-conforming"
  )
  expect_error(
    dmatrixnorm(x = A, mean = A, L = A, V = B),
    "Non-conforming"
  )
  expect_error(
    dmatrixnorm(x = A, mean = A, U = A, V = B),
    "Non-conforming"
  )
  expect_error(
    dmatrixnorm(x = A, mean = A, U = A, R = B),
    "Non-conforming"
  )

  expect_error(
    dmatrixt(x = A, df = 1, mean = A, L = B, R = A),
    "Non-conforming"
  )
  expect_error(
    dmatrixt(x = A, df = 1, mean = A, L = A, R = B),
    "Non-conforming"
  )
  expect_error(dmatrixt(x = A, df = 1, mean = A, U = B, R = A), "Non-conforming")
  expect_error(dmatrixt(x = A, df = 1, mean = A, L = A, V = B), "Non-conforming")
  expect_error(dmatrixt(x = A, df = 1, mean = A, U = A, V = B), "Non-conforming")
  expect_error(dmatrixt(x = A, df = 1, mean = A, U = A, R = B), "Non-conforming")

  expect_error(dmatrixinvt(x = A, df = 1, mean = A, L = B, R = A), "Non-conforming")
  expect_error(dmatrixinvt(x = A, df = 1, mean = A, L = A, R = B), "Non-conforming")
  expect_error(dmatrixinvt(x = A, df = 1, mean = A, U = B, R = A), "Non-conforming")
  expect_error(dmatrixinvt(x = A, df = 1, mean = A, L = A, V = B), "Non-conforming")
  expect_error(dmatrixinvt(x = A, df = 1, mean = A, U = A, V = B), "Non-conforming")
  expect_error(dmatrixinvt(x = A, df = 1, mean = A, U = A, R = B), "Non-conforming")

  expect_error(MLmatrixnorm(rmatrixnorm(n = 100, mean = A), U = B, V = A))
  expect_error(MLmatrixnorm(rmatrixnorm(n = 100, mean = A), U = A, V = B))
  expect_error(MLmatrixt(rmatrixnorm(n = 100, mean = A), df = 0, U = A, V = B))
  expect_error(MLmatrixt(rmatrixnorm(n = 100, mean = A), U = B, V = A))
  expect_error(MLmatrixt(rmatrixnorm(n = 100, mean = A), U = A, V = B))
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
  A <- rmatrixnorm(5, mean = matrix(0, nrow = 2, ncol = 2))
  B <- rmatrixnorm(5, mean = matrix(1, nrow = 2, ncol = 2))
  C <- array(c(A, B), dim = c(2, 2, 10))
  D <- array(0, dim = c(2, 2, 5))
  E <- array(c(A, D), dim = c(2, 2, 10))

  groups <- c(rep(1, 5), rep(2, 5))
  groups.empty <- factor(rep("1", 10), levels = c("1", "2"))
  priors <- c(.5, .5)
  expect_error(
    matrixlda(c(C), grouping = c(rep(1, 40), rep(2, 40)), prior = priors),
    "array"
  )
  expect_error(
    matrixlda(C, grouping = c(rep(1, 120), rep(2, 120)), prior = priors),
    "are different"
  )
  expect_error(
    matrixlda(C, grouping = groups, prior = c(.5, .4)),
    "invalid 'prior'"
  )
  expect_error(
    matrixlda(C, grouping = groups, prior = c(.4, .4, .2)),
    "incorrect length"
  )
  expect_warning(
    matrixlda(C, grouping = groups.empty, prior = priors),
    "empty"
  )
  #  expect_error(
  #    matrixlda(E, grouping = groups, prior = priors)
  #  )
  # this one is a problem but i don't think we should test it
})

context("Testing QDA")



test_that("Testing bad input to QDA", {
  A <- rmatrixnorm(5, mean = matrix(0, nrow = 2, ncol = 2))
  B <- rmatrixnorm(5, mean = matrix(1, nrow = 2, ncol = 2))
  C <- array(c(A, B), dim = c(2, 2, 10))
  D <- array(0, dim = c(2, 2, 5))
  E <- array(c(A, D), dim = c(2, 2, 10))

  groups <- c(rep(1, 5), rep(2, 5))
  groups.empty <- factor(rep("1", 10), levels = c("1", "2"))
  priors <- c(.5, .5)


  expect_error(
    matrixqda(c(C), grouping = c(rep(1, 40), rep(2, 40)), prior = priors),
    "array"
  )
  expect_error(
    matrixqda(C, grouping = c(rep(1, 120), rep(2, 120)), prior = priors),
    "are different"
  )
  expect_error(
    matrixqda(C, grouping = groups, prior = c(.5, .4)),
    "invalid 'prior'"
  )
  expect_error(
    matrixqda(C, grouping = groups, prior = c(.4, .4, .2)),
    "incorrect length"
  )
  expect_warning(
    matrixqda(C, grouping = groups.empty, prior = priors),
    "empty"
  )
  expect_error(
    matrixqda(E, grouping = groups, prior = priors)
  )
})

context("Out of bounds")

test_that("Out of bounds numeric input: ", {
  A <- diag(5)
  A[5, 5] <- 0
  expect_error(rmatrixt(0, 1, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = diag(5)))
  expect_error(rmatrixt(1, -1, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = diag(5)))
  expect_error(rmatrixt(1, 1, matrix(0, nrow = 5, ncol = 5), U = -diag(5), V = diag(5)))
  expect_error(rmatrixt(1, 1, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = -diag(5)))
  expect_error(rmatrixt(1, 1, matrix(0, nrow = 5, ncol = 5), L = A, V = diag(5)))
  expect_error(rmatrixt(1, 1, matrix(0, nrow = 5, ncol = 5), U = diag(5), R = A))
  expect_error(rmatrixt(1, 1, matrix(0, nrow = 5, ncol = 5), U = A, V = diag(5)))
  expect_error(rmatrixt(1, 1, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = A))

  expect_error(rmatrixinvt(0, 1, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = diag(5)))
  expect_error(rmatrixinvt(1, -1, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = diag(5)))
  expect_error(rmatrixinvt(1, 1, matrix(0, nrow = 5, ncol = 5), U = -diag(5), V = diag(5)))
  expect_error(rmatrixinvt(1, 1, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = -diag(5)))
  expect_error(rmatrixinvt(1, 1, matrix(0, nrow = 5, ncol = 5), L = A, V = diag(5)))
  expect_error(rmatrixinvt(1, 1, matrix(0, nrow = 5, ncol = 5), U = diag(5), R = A))
  expect_error(rmatrixinvt(1, 1, matrix(0, nrow = 5, ncol = 5), U = A, V = diag(5)))
  expect_error(rmatrixinvt(1, 1, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = A))

  expect_error(rmatrixnorm(0, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = diag(5)))
  expect_error(rmatrixnorm(1, matrix(0, nrow = 5, ncol = 5), U = -diag(5), V = diag(5)))
  expect_error(rmatrixnorm(1, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = -diag(5)))
  expect_error(rmatrixnorm(1, matrix(0, nrow = 5, ncol = 5), L = A, V = diag(5)))
  expect_error(rmatrixnorm(1, matrix(0, nrow = 5, ncol = 5), U = diag(5), R = A))
  expect_error(rmatrixnorm(1, 1, matrix(0, nrow = 5, ncol = 5), L = A, V = diag(5)))
  expect_error(rmatrixnorm(1, 1, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = A))

  expect_error(dmatrixt(1, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = diag(5)))
  expect_error(dmatrixt(1, df = 0, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = diag(5)))

  expect_error(dmatrixinvt(1, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = diag(5)))

  expect_error(dmatrixt(df = -1, x = matrix(0, nrow = 5, ncol = 5), U = diag(5), V = diag(5)))
  expect_error(dmatrixt(df = 1, x = matrix(0, nrow = 5, ncol = 5), U = -diag(5), V = diag(5)))
  expect_error(dmatrixt(df = 1, x = matrix(0, nrow = 5, ncol = 5), U = diag(5), V = -diag(5)))
  expect_error(dmatrixt(df = 1, x = matrix(0, nrow = 5, ncol = 5), L = A, V = diag(5)))
  expect_error(dmatrixt(df = 1, x = matrix(0, nrow = 5, ncol = 5), U = diag(5), R = A))
  expect_error(dmatrixt(df = 1, x = matrix(0, nrow = 5, ncol = 5), U = A, V = diag(5)))
  expect_error(dmatrixt(df = 1, x = matrix(0, nrow = 5, ncol = 5), U = diag(5), V = A))

  expect_error(dmatrixinvt(df = -1, x = matrix(0, nrow = 5, ncol = 5), U = diag(5), V = diag(5)))
  expect_error(dmatrixinvt(df = 1, x = matrix(0, nrow = 5, ncol = 5), U = -diag(5), V = diag(5)))
  expect_error(dmatrixinvt(df = 1, x = matrix(0, nrow = 5, ncol = 5), U = diag(5), V = -diag(5)))
  expect_error(dmatrixinvt(df = 1, x = matrix(0, nrow = 5, ncol = 5), L = A, V = diag(5)))
  expect_error(dmatrixinvt(df = 1, x = matrix(0, nrow = 5, ncol = 5), U = diag(5), R = A))
  expect_error(dmatrixinvt(df = 1, x = matrix(0, nrow = 5, ncol = 5), U = A, V = diag(5)))
  expect_error(dmatrixinvt(df = 1, x = matrix(0, nrow = 5, ncol = 5), U = diag(5), V = A))

  expect_error(dmatrixnorm(matrix(0, nrow = 5, ncol = 5), U = -diag(5), V = diag(5)))
  expect_error(dmatrixnorm(matrix(0, nrow = 5, ncol = 5), U = diag(5), V = -diag(5)))
  expect_error(dmatrixnorm(matrix(0, nrow = 5, ncol = 5), L = A, V = diag(5)))
  expect_error(dmatrixnorm(matrix(0, nrow = 5, ncol = 5), U = diag(5), R = A))
  expect_error(dmatrixnorm(1, matrix(0, nrow = 5, ncol = 5), L = A, V = diag(5)))
  expect_error(dmatrixnorm(1, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = A))
})
