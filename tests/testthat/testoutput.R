library(matrixdist)

context("Checking outputs match")


test_that("Testing helper functions:", {

  # A = diag(5) + 1
  # B = posmatsqrt(A)
  # C = posmatsqrtinv(A)
  #
  # expect_equal(B %*% C, diag(5))
  # expect_equal(B, t(B))
  # expect_equal(C, t(C))
  # expect_equal(A, (B %*% B))
  C = matrix(c(1, .5, .25, .5, 1, .5, .25, .5, 1), nrow = 3)
  expect_equal(ARgenerate(3, .5), C)


})

test_that("Equivalent outputs for different options:", {
  set.seed(2018020201)
  A <-
    rmatrixnorm(
      n = 1,
      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
      L = matrix(c(2, 1, 0, .1), nrow = 2),
      list = FALSE
    )
  set.seed(2018020201)
  B <-
    rmatrixnorm(
      n = 10,
      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
      L = matrix(c(2, 1, 0, .1), nrow = 2),
      list = TRUE
    )
  set.seed(2018020201)
  C <-
    rmatrixnorm(
      n = 10,
      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
      L = matrix(c(2, 1, 0, .1), nrow = 2),
      list = FALSE
    )
  expect_equal(A, B[[1]])
  expect_equal(A, C[, , 1])
  set.seed(2018020202)
  A <-
    rmatrixt(
      n = 1,
      df = 2,
      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
      L = matrix(c(2, 1, 0, .1), nrow = 2),
      list = FALSE
    )
  set.seed(2018020202)
  B <-
    rmatrixt(
      n = 1,
      df = 2,
      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
      L = matrix(c(2, 1, 0, .1), nrow = 2),
      list = TRUE
    )
  set.seed(2018020202)
  C <-
    rmatrixt(
      n = 1,
      df = 2,
      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
      L = matrix(c(2, 1, 0, .1), nrow = 2),
      array = TRUE
    )
  expect_equal(A, B[[1]])
  expect_equal(A, C[, , 1])
  set.seed(2018020203)
  A <-
    rmatrixinvt(
      n = 1,
      df = 2,
      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
      L = matrix(c(2, 1, 0, .1), nrow = 2),
      list = FALSE
    )
  set.seed(2018020203)
  B <-
    rmatrixinvt(
      n = 1,
      df = 2,
      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
      L = matrix(c(2, 1, 0, .1), nrow = 2),
      list = TRUE
    )
  set.seed(2018020203)
  C <-
    rmatrixinvt(
      n = 1,
      df = 2,
      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
      L = matrix(c(2, 1, 0, .1), nrow = 2),
      array = TRUE
    )
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

  expect_equal(dnorm(1), dmatrixnorm(matrix(1)))
  expect_equal(dnorm(1), dmatrixnorm.unroll(matrix(1)))

  expect_equal(dmatrixnorm.unroll(x = matrix(.5, nrow = 10, ncol = 2),
                                  log = TRUE),
               dmatrixnorm(x = matrix(.5, nrow = 10,ncol = 2),
                           log = TRUE))
  expect_equal(dmatrixnorm.unroll(x = matrix(.5, nrow = 2, ncol = 10),
                                  log = TRUE),
               dmatrixnorm(x = matrix(.5, nrow = 10,ncol = 2),
                           log = TRUE))


  expect_equal(dt(1, 1), (dmatrixt(1, 1)))
  expect_equal(dmatrixt(matrix(1), df = 10, U = 10 * matrix(1)),
               dt(1, 10),
               tolerance = 1e-6)
  expect_equal(dmatrixt(matrix(1), df = 10, V = 10 * matrix(1)),
               dt(1, 10),
               tolerance = 1e-6)
  A = as.vector(rmatrixnorm(1e4, 0, list = FALSE))
  B = as.vector(rnorm(1e4))
  expect_equal(var(A), var(B), tolerance = .08)
  expect_equal(mean(A), mean(B), tolerance = .035)
  df = 10
  A = as.vector(rmatrixt(
    1e4,
    df = df,
    0,
    V = df,
    list = FALSE
  ))
  B = as.vector(rt(1e4, df = 10))
  expect_equal(var(A), var(B), tolerance = .08)
  expect_equal(mean(A), mean(B), tolerance = .035)
  U.one = V.two = matrix(1)
  dim.one = c(1, 6)
  dim.two = c(6, 1)
  x = array(rep(1, 6), dim = c(1, 6))
  U.two = V.one = ARgenerate(6, .7)
  df = 5

  #dmvt(x,sigma = U.two,df = 5)
  expect_equal(
    dmatrixt(
      x,
      df,
      U = U.one,
      V = df * V.one,
      log = T
    ),
    dmatrixt(
      t(x),
      df,
      U = U.two,
      V = df * V.two,
      log = T
    ),
    tolerance = .000001
  )
  expect_equal(dmatrixt(x, df, U = U.one, V = df * V.one, log = T),
  -4.663386,
  tolerance = .000001)
  expect_equal( dmatrixt(t(rep(1,5)),df = 5,U = 5,log = TRUE),
                dmatrixt((rep(1,5)),df = 5,V = 5,log = TRUE))
  expect_equal(dmatrixt((rep(1,5)),df = 5,V = 5,log = TRUE),
               -7.457784, tolerance = 1e-6)

  expect_equal( dmatrixinvt(t(rep(.5,5)),df = 10,U = 5,log = TRUE),
                 dmatrixinvt((rep(.5,5)),df = 10,V = 5,log = TRUE))

  set.seed(20180222)
  A <- rWishart(1,7,diag(6))[,,1]
  expect_equal(dmatrixt(t(rep(1,6)), df = 5, U = 5, V = A, log = TRUE),
               dmatrixt((rep(1,6)),df = 5,V = 5,U = A, log = TRUE))
  expect_equal(dmatrixt(t(rep(1,6)), df = 5, U = 5, V = A, log = TRUE),
          -16.07342, tolerance = 1e-6)


  set.seed(20180219)
  A <- rmatrixnorm(40, mean = array(1, dim = c(4, 5)),
                U = CSgenerate(4,.2), V = ARgenerate(5, .8))
  B <- MLmatrixnorm(A, row.mean = T)
  expect_equal(B$U[1, 1], 1)
  expect_equal(B$V[1, 1], 1)
  expect_true(B$convergence)
  expect_equal(B$mean[1, 1], B$mean[1, 2])
  C <- MLmatrixnorm(A, col.mean = T)
  expect_equal(C$mean[1, 1], C$mean[2, 1])
  C <- MLmatrixnorm(A, col.variance = "CS")
  expect_equal(C$V[1,2],C$V[1,5])
  C <- MLmatrixnorm(A, row.variance = "AR(1)")
  expect_equal(C$U[2,1],C$U[3,2])
  C <- MLmatrixnorm(A, row.variance = "I")
  expect_equal(C$U[1,2],0)
  C <- MLmatrixnorm(A, col.variance = "I")
  expect_equal(C$V[1,4],0)
})


test_that("Output of LDA/QDA/Predict", {
  A <- rmatrixnorm(4, mean = matrix(0, nrow = 2, ncol = 2))
  B <- rmatrixnorm(4, mean = matrix(1, nrow = 2, ncol = 2))
  C <- array(c(A,B), dim = c(2,2,8))
  D <- array(0, dim = c(2,2,4))
  E <- array(c(A,D), dim = c(2,2,8))
  groups <- c(rep(1,4),rep(2,4))
  groups.empty <- factor(rep("1",8), levels = c("1","2"))
  priors = c(.5,.5)
  ldamodel <- matrixlda(C, groups, priors)
  qdamodel <- matrixqda(C, groups, priors)
  expect_error(predict(ldamodel, newdata = matrix(0,nrow = 3, ncol = 2)),
               "dimension")
  expect_error(predict(ldamodel, newdata = matrix(0,nrow = 2, ncol = 3)),
               "dimension")
  expect_error(predict(qdamodel, newdata = matrix(0,nrow = 3, ncol = 2)),
               "dimension")
  expect_error(predict(qdamodel, newdata = matrix(0,nrow = 2, ncol = 3)),
               "dimension")


  expect_equal(sum(predict(ldamodel, newdata = matrix(
    0, nrow = 2, ncol = 2))$posterior), 1)

  expect_equal(sum(predict(qdamodel, newdata = matrix(
    0, nrow = 2, ncol = 2))$posterior), 1)

  newlda <- ldamodel
  newqda <- qdamodel
  newprior <- c(-1,2)

  expect_error(predict(newlda, prior = newprior),
               "invalid 'prior'")
  expect_error(predict(newqda, prior = newprior),
               "invalid 'prior'")
  newprior <- c(.7,.7)
  expect_error(predict(newlda, prior = newprior),
               "invalid 'prior'")
  expect_error(predict(newqda, prior = newprior),
               "invalid 'prior'")


})
