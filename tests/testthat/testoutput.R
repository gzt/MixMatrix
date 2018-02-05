library(matrixdist)

context("Checking outputs match")


test_that("Testing helper functions:",{
  expect_equal(lgamma(1:12)-lmvgamma(1:12,1),array(0,dim=12),tolerance=1e-6)
  expect_equal(gamma(1:12) - mvgamma(1:12,1),array(0,dim=12),tolerance=1e-6)
  A = diag(5) + 1
  B = posmatsqrt(A)
  expect_equal(B,t(B))
  expect_equal(A, (B %*% B))
})
