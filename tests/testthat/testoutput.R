library(matrixdist)

context("Checking outputs match")


test_that("Testing helper functions:", {
  expect_equal(lgamma(1:12) - lmvgamma(1:12, 1), array(0, dim = 12), tolerance = 1e-7)
  expect_equal(gamma(1:12) - mvgamma(1:12, 1), array(0, dim = 12), tolerance = 1e-7)
  p = 2
  expect_equal( (p * (p - 1)/4 * log(pi) + lgamma(5 - 0) + lgamma(5 - .5)),as.numeric(lmvgamma(5,2)))
  A = diag(5) + 1
  B = posmatsqrt(A)
  C = posmatsqrtinv(A)

  expect_equal(B %*% C, diag(5))
  expect_equal(B, t(B))
  expect_equal(C, t(C))
  expect_equal(A, (B %*% B))
  C = matrix(c(1, .5, .25, .5, 1, .5, .25, .5, 1), nrow = 3)
  expect_equal(toepgenerate(3, .5), C)
})

test_that("Equivalent outputs for different options:", {
  set.seed(2018020201)
  A <- rmatrixnorm(n = 1, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
     L = matrix(c(2, 1, 0, .1), nrow = 2), list = FALSE)
  set.seed(2018020201)
  B <- rmatrixnorm(n = 10, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
     L = matrix(c(2, 1, 0, .1), nrow = 2), list = TRUE)
  set.seed(2018020201)
  C <- rmatrixnorm(n = 10, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
     L = matrix(c(2, 1, 0, .1), nrow = 2), list = FALSE)
  expect_equal(A, B[[1]])
  expect_equal(A, C[, , 1])
  set.seed(2018020202)
  A <- rmatrixt(n = 1, df = 2, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
                   L = matrix(c(2, 1, 0, .1), nrow = 2), list = FALSE)
  set.seed(2018020202)
  B <- rmatrixt(n = 1, df = 2, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
                   L = matrix(c(2, 1, 0, .1), nrow = 2), list = TRUE)
  set.seed(2018020202)
  C <- rmatrixt(n = 1, df = 2, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
                   L = matrix(c(2, 1, 0, .1), nrow = 2), array = TRUE)
  expect_equal(A, B[[1]])
  expect_equal(A, C[, , 1])
  set.seed(2018020203)
  A <- rmatrixinvt(n = 1, df = 2, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
                L = matrix(c(2, 1, 0, .1), nrow = 2), list = FALSE)
  set.seed(2018020203)
  B <- rmatrixinvt(n = 1, df = 2, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
                L = matrix(c(2, 1, 0, .1), nrow = 2), list = TRUE)
  set.seed(2018020203)
  C <- rmatrixinvt(n = 1, df = 2, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
                L = matrix(c(2, 1, 0, .1), nrow = 2), array = TRUE)
  expect_equal(A, B[[1]])
  expect_equal(A, C[, , 1])

})

test_that("Equivalent outputs for different functions:", {
set.seed(20180205)
  A = as.vector(rmatrixnorm(10, 0, array = TRUE))
set.seed(20180205)
  B = rnorm(10)
set.seed(20180205)
  C = as.vector(rmatrixt(10, df = 0, 0, array = TRUE))
set.seed(20180205)
  D = as.vector(rmatrixt(10, df = Inf, 0, array = TRUE))
  expect_equal(A, B)
  expect_equal(A, C)
  expect_equal(A, D)

  expect_equal(dnorm(1), dmatrixnorm(1))
  expect_equal(dt(1, 1), (dmatrixt(1, 1)))

  A = as.vector(rmatrixnorm(1e4, 0, list = FALSE))
  B = as.vector(rnorm(1e4))
  expect_equal(var(A), var(B), tolerance = .08)
  expect_equal(mean(A), mean(B), tolerance = .035)

  A = as.vector(rmatrixt(1e4,df = 10, 0, list = FALSE))
  B = as.vector(rt(1e4, df = 10))
  expect_equal(var(A), var(B), tolerance = .08)
  expect_equal(mean(A), mean(B), tolerance = .035)
  U.one = V.two = matrix(1)
  dim.one = c(1,6)
  dim.two = c(6,1)
  x = array(rep(1,6),dim = c(1,6))
  U.two = V.one = toepgenerate(6, .7)
  df = 5

  #dmvt(x,sigma = U.two,df = 5)
  expect_equal(dmatrixt(x, df, U = U.one, V = V.one, log = T),
               dmatrixt(t(x), df, U = U.two, V = V.two, log = T),
               tolerance = .0001)
  expect_equal(dmatrixt(x, df, U = U.one, V = V.one, log = T),
               -4.663386, tolerance = .0001)
  set.seed(20180211)
  A <- rInvCholWishart(1,10,diag(5))[,,1]
  set.seed(20180211)
  B <- rCholWishart(1,10,diag(5))[,,1]
  set.seed(20180211)
  C <- chol(rWishart(1,10,diag(5))[,,1])

  expect_equal(sum(abs(A[lower.tri(A)])),0)
  expect_equal(sum(abs(B[lower.tri(B)])),0)
  expect_equal(A %*% B, diag(5))
  expect_equal(B,C)

})
