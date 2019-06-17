library(MixMatrix)

context("Testing input positive definite/invertible integrity")


test_that("Bad rank in covariance:", {
  A = matrix(0, nrow = 2, ncol = 2)

  expect_error(rmatrixnorm(
    2,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    L = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(rmatrixnorm(
    2,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    R = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_length(rmatrixnorm(
    2,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    L = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2),
    force = TRUE
  ), 8)
  expect_length(rmatrixnorm(
    2,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    R = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2),
    force = TRUE
  ), 8)
  expect_error(rmatrixnorm(
    2,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    U = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(rmatrixnorm(
    2,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    V = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))

  expect_error(rmatrixt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    L = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(rmatrixt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    R = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_length(rmatrixt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    R = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2),
    force = TRUE
  ), 8)
  expect_error(rmatrixt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    U = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(rmatrixt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    V = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))


  expect_error(rmatrixinvt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    L = matrix(c(1, 1, .5, .5), nrow = 2, ncol = 2)
  ))
  expect_error(rmatrixinvt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    R = matrix(c(1, 1, .5, .5), nrow = 2, ncol = 2)
  ))
  expect_error(rmatrixinvt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    U = matrix(c(1, 1, .5, .5), nrow = 2, ncol = 2)
  ))
  expect_error(rmatrixinvt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    V = matrix(c(1, 1, .5, .5), nrow = 2, ncol = 2)
  ))

  expect_error(dmatrixnorm(
    A,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    L = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixnorm(
    A,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    R = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixnorm(
    A,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    U = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixnorm(
    A,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    V = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))

  expect_error(dmatrixinvt(
    A,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    L = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixinvt(
    A,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    R = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixinvt(
    A,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    U = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixinvt(
    A,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    V = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))

  expect_error(dmatrixt(
    A,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    L = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixt(
    A,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    R = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixt(
    A,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    U = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixt(
    A,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    V = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))



})
