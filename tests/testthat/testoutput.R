library(MixMatrix)

context("Checking outputs match")


test_that("Testing helper functions:", {
  c_mat <- matrix(c(1, .5, .25, .5, 1, .5, .25, .5, 1), nrow = 3)
  expect_equal(ARgenerate(3, .5), c_mat)
})

test_that("Equivalent outputs for different options:", {
  set.seed(2018020201)
  a_mat <- rmatrixnorm(
    n = 1, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c(2, 1, 0, .1), nrow = 2),
    list = FALSE
  )
  set.seed(2018020201)
  b_mat <- rmatrixnorm(
    n = 10, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c(2, 1, 0, .1), nrow = 2), list = TRUE
  )
  set.seed(2018020201)
  c_mat <- rmatrixnorm(
    n = 10, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c(2, 1, 0, .1), nrow = 2),
    list = FALSE
  )
  expect_equal(a_mat, b_mat[[1]])
  expect_equal(a_mat, c_mat[, , 1])
  expect_equal(
    dmatrixnorm(a_mat, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2)),
    dmatrixnorm(b_mat[[1]], mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2))
  )
  expect_equal(
    dmatrixnorm(a_mat, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2)),
    dmatrixnorm(c_mat, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2))[1]
  )

  set.seed(2018020202)
  a_mat <- rmatrixt(
    n = 1, df = 2, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c(2, 1, 0, .1), nrow = 2),
    list = FALSE
  )
  set.seed(2018020202)
  b_mat <- rmatrixt(
    n = 1, df = 2, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c(2, 1, 0, .1), nrow = 2),
    list = TRUE
  )
  set.seed(2018020202)
  c_mat <- rmatrixt(
    n = 1, df = 2, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c(2, 1, 0, .1), nrow = 2),
    array = TRUE
  )

  expect_equal(a_mat, b_mat[[1]])
  expect_equal(a_mat, c_mat[, , 1])
  expect_equal(
    dmatrixt(a_mat, df = 2, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2)),
    dmatrixt(b_mat[[1]],
      df = 2,
      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2)
    )
  )
  expect_equal(
    dmatrixt(a_mat, df = 2, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2)),
    dmatrixt(c_mat, df = 2, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2))
  )

  set.seed(2018020203)
  a_mat <- rmatrixinvt(
    n = 1, df = 2, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c(2, 1, 0, .1), nrow = 2),
    list = FALSE
  )
  set.seed(2018020203)
  b_mat <- rmatrixinvt(
    n = 1, df = 2, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c(2, 1, 0, .1), nrow = 2),
    list = TRUE
  )
  set.seed(2018020203)
  c_mat <- rmatrixinvt(
    n = 1, df = 2,
    mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c(2, 1, 0, .1), nrow = 2),
    array = TRUE
  )
  expect_equal(a_mat, b_mat[[1]])
  expect_equal(a_mat, c_mat[, , 1])

  expect_equal(
    dmatrixnorm(array(1:10, dim = c(1, 1, 10)), 0, U = 4),
    dnorm(1:10, 0, 2)
  )
  expect_equal(
    dmatrixnorm(array(1:10, dim = c(1, 1, 10)), 0, V = 16),
    dnorm(1:10, 0, 4)
  )
})

test_that("Equivalent outputs for different functions:", {
  set.seed(20180205)
  a_mat <- as.vector(rmatrixnorm(10, 0, array = TRUE))
  set.seed(20180205)
  b_mat <- rnorm(10)
  set.seed(20180205)
  c_mat <- as.vector(rmatrixt(10, df = 0, 0, array = TRUE))
  set.seed(20180205)
  d_mat <- as.vector(rmatrixt(10, df = Inf, 0, array = TRUE))
  expect_equal(a_mat, b_mat)
  expect_equal(a_mat, c_mat)
  expect_equal(a_mat, d_mat)

  expect_equal(dnorm(1), dmatrixnorm(matrix(1)))
  expect_equal(dnorm(1), dmatrixnorm.unroll(matrix(1)))

  expect_equal(
    dmatrixnorm.unroll(
      x = matrix(.5, nrow = 10, ncol = 2),
      log = TRUE
    ),
    dmatrixnorm(
      x = matrix(.5, nrow = 10, ncol = 2),
      log = TRUE
    )
  )
  expect_equal(
    dmatrixnorm.unroll(
      x = matrix(.5, nrow = 2, ncol = 10),
      log = TRUE
    ),
    dmatrixnorm(
      x = matrix(.5, nrow = 10, ncol = 2),
      log = TRUE
    )
  )


  expect_equal(dt(1, 1), (dmatrixt(1, 1)))
  expect_equal(dmatrixt(matrix(1), df = 10, U = 10 * matrix(1)),
    dt(1, 10),
    tolerance = 1e-6
  )
  expect_equal(dmatrixt(matrix(1), df = 10, V = 10 * matrix(1)),
    dt(1, 10),
    tolerance = 1e-6
  )
  a_mat <- as.vector(rmatrixnorm(1e4, 0, list = FALSE))
  b_mat <- as.vector(rnorm(1e4))
  expect_equal(var(a_mat), var(b_mat), tolerance = .08)
  expect_equal(mean(a_mat), mean(b_mat), tolerance = .035)
  df <- 10
  a_mat <- as.vector(rmatrixt(1e4, df = df, 0, V = df, list = FALSE))
  b_mat <- as.vector(rt(1e4, df = 10))
  expect_equal(var(a_mat), var(b_mat), tolerance = .08)
  expect_equal(mean(a_mat), mean(b_mat), tolerance = .035)
  u_one <- v_two <- matrix(1)
  dim.one <- c(1, 6)
  dim.two <- c(6, 1)
  x <- array(rep(1, 6), dim = c(1, 6))
  u_two <- v_one <- ARgenerate(6, .7)
  df <- 5

  expect_equal(dmatrixt(x, df, U = u_one, V = df * v_one, log = T),
    dmatrixt(t(x), df, U = u_two, V = df * v_two, log = T),
    tolerance = .000001
  )
  expect_equal(dmatrixt(x, df, U = u_one, V = df * v_one, log = T),
    -4.663386,
    tolerance = .000001
  )
  expect_equal(
    dmatrixt(t(rep(1, 5)), df = 5, U = 5, log = TRUE),
    dmatrixt((rep(1, 5)), df = 5, V = 5, log = TRUE)
  )
  expect_equal(dmatrixt((rep(1, 5)), df = 5, V = 5, log = TRUE),
    -7.457784,
    tolerance = 1e-6
  )

  expect_equal(
    dmatrixinvt(t(rep(.5, 5)), df = 10, U = 5, log = TRUE),
    dmatrixinvt((rep(.5, 5)), df = 10, V = 5, log = TRUE)
  )

  set.seed(20180222)
  a_mat <- rWishart(1, 7, diag(6))[, , 1]
  expect_equal(
    dmatrixt(t(rep(1, 6)), df = 5, U = 5, V = a_mat, log = TRUE),
    dmatrixt((rep(1, 6)), df = 5, V = 5, U = a_mat, log = TRUE)
  )
  expect_equal(dmatrixt(t(rep(1, 6)), df = 5, U = 5, V = a_mat, log = TRUE),
    -16.07342,
    tolerance = 1e-6
  )


  set.seed(20180219)
  a_mat <- rmatrixnorm(40,
    mean = array(1, dim = c(4, 5)),
    U = CSgenerate(4, .2), V = ARgenerate(5, .8), list = TRUE
  )
  b_mat <- MLmatrixnorm(a_mat, row.mean = TRUE, V = ARgenerate(5, .8))
  expect_equal(b_mat$U[1, 1], 1)
  expect_equal(b_mat$V[1, 1], 1)

  expect_true(b_mat$convergence)
  expect_equal(b_mat$mean[1, 1], b_mat$mean[1, 2])
  c_mat <- MLmatrixnorm(a_mat,
    col.mean = TRUE,
    row.mean = TRUE, U = CSgenerate(4, .2)
  )
  expect_equal(c_mat$mean[1, 1], c_mat$mean[2, 1])
  expect_equal(c_mat$mean[1, 1], c_mat$mean[1, 2])


  c_mat <- MLmatrixnorm(a_mat, col.mean = TRUE)
  expect_equal(c_mat$mean[1, 1], c_mat$mean[2, 1])
  c_mat <- MLmatrixnorm(a_mat, row.variance = "CS")
  expect_equal(c_mat$U[1, 2], c_mat$U[1, 4])
  c_mat <- MLmatrixnorm(a_mat, col.variance = "CS")
  expect_equal(c_mat$V[1, 2], c_mat$V[1, 4])
  c_mat <- MLmatrixnorm(a_mat, row.variance = "AR(1)")
  expect_equal(c_mat$U[2, 1], c_mat$U[3, 2])
  c_mat <- MLmatrixnorm(a_mat, col.variance = "AR(1)")
  expect_equal(c_mat$V[2, 1], c_mat$V[3, 2])
  c_mat <- MLmatrixnorm(a_mat, row.variance = "corr")
  expect_equal(c_mat$U[2, 2], c_mat$U[1, 1])
  c_mat <- MLmatrixnorm(a_mat, col.variance = "corr")
  expect_equal(c_mat$V[2, 2], c_mat$V[3, 3])
  c_mat <- MLmatrixnorm(a_mat, row.variance = "I")
  expect_equal(c_mat$U[1, 2], 0)
  c_mat <- MLmatrixnorm(a_mat, col.variance = "I")
  expect_equal(c_mat$V[1, 4], 0)

  d_mat <- MLmatrixt(a_mat,
    col.mean = TRUE,
    U = CSgenerate(4, .2), V = ARgenerate(5, .8)
  )
  expect_true(d_mat$convergence)
  expect_warning(MLmatrixt(a_mat, fixed = FALSE, max.iter = 2))
  expect_equal(d_mat$U[1, 1], 1)
  expect_equal(d_mat$V[1, 1], 1)
  expect_equal(d_mat$mean[1, 1], d_mat$mean[2, 1])
  c_mat <- MLmatrixt(a_mat, col.mean = TRUE, row.mean = TRUE)
  expect_equal(c_mat$mean[1, 1], c_mat$mean[2, 1])
  expect_equal(c_mat$mean[1, 1], c_mat$mean[1, 2])

  c_mat <- MLmatrixt(a_mat, col.variance = "CS")
  expect_equal(c_mat$V[1, 2], c_mat$V[1, 5])
  c_mat <- MLmatrixt(a_mat, row.variance = "CS")
  expect_equal(c_mat$U[1, 2], c_mat$U[1, 4])
  c_mat <- MLmatrixt(a_mat, row.variance = "AR(1)")
  expect_equal(c_mat$U[2, 1], c_mat$U[3, 2])
  c_mat <- MLmatrixt(a_mat, col.variance = "AR(1)")
  expect_equal(c_mat$V[2, 1], c_mat$V[3, 2])
  c_mat <- MLmatrixt(a_mat, col.variance = "corr")
  expect_equal(c_mat$V[1, 1], c_mat$V[2, 2])
  c_mat <- MLmatrixt(a_mat, row.variance = "corr")
  expect_equal(c_mat$U[1, 1], c_mat$U[2, 2])
  c_mat <- MLmatrixt(a_mat, row.variance = "I")
  expect_equal(c_mat$U[1, 2], 0)
  c_mat <- MLmatrixt(a_mat, col.variance = "I")
  expect_equal(c_mat$V[1, 4], 0)
})

context("Testing LDA/QDa_mat output")
test_that("Output of LDA/QDA/Predict", {
  set.seed(20190628)
  a_mat <- rmatrixnorm(4, mean = matrix(0, nrow = 2, ncol = 2))
  b_mat <- rmatrixnorm(4, mean = matrix(1, nrow = 2, ncol = 2))
  set.seed(20190628)

  c_mat <- array(c(a_mat, b_mat), dim = c(2, 2, 8))
  c_matzero <- c_mat
  c_matzero[1, 1, ] <- 0
  d_mat <- array(0, dim = c(2, 2, 4))
  e_mat <- array(c(a_mat, d_mat), dim = c(2, 2, 8))
  groups <- c(rep(1, 4), rep(2, 4))

  priors <- c(.5, .5)
  ldamodel <- matrixlda(c_mat, groups, priors, subset = rep(TRUE, 8))
  qdamodel <- matrixqda(c_mat, groups, priors, subset = rep(TRUE, 8))
  expect_error(matrixlda(c_matzero, groups, priors,
    subset = rep(TRUE, 8)
  ), "constant")
  expect_error(suppressWarnings(matrixqda(c_matzero, groups, priors,
    subset = rep(TRUE, 8)
  )), "constant")
  expect_error(
    predict(ldamodel, newdata = matrix(0, nrow = 3, ncol = 2)),
    "dimension"
  )
  expect_error(
    predict(ldamodel, newdata = (matrix(0, nrow = 2, ncol = 3))),
    "dimension"
  )
  expect_error(
    predict(qdamodel, newdata = matrix(0, nrow = 3, ncol = 2)),
    "dimension"
  )
  expect_error(
    predict(qdamodel, newdata = (matrix(0, nrow = 2, ncol = 3))),
    "dimension"
  )

  expect_equal(sum(predict(ldamodel, newdata = matrix(
    0,
    nrow = 2, ncol = 2
  ))$posterior), 1)
  expect_equal(sum(predict(ldamodel, prior = c(.7, .3))$posterior[1, ]), 1)

  expect_equal(sum(predict(qdamodel, prior = c(.7, .3), newdata = matrix(
    0,
    nrow = 2, ncol = 2
  ))$posterior), 1)
  expect_equal(sum(predict(qdamodel)$posterior[1, ]), 1)

  llda <- logLik(ldamodel)
  lqda <- logLik(qdamodel)
  expect_equal(class(llda), "logLik")
  expect_equal(nobs(llda), 8)
  expect_equal(class(lqda), "logLik")
  expect_equal(nobs(lqda), 8)
  expect_equal(attributes(llda)$nobs, nobs(ldamodel))
  expect_equal(attributes(lqda)$nobs, nobs(qdamodel))
  newlda <- matrixlda(c_mat, groups, priors, method = "t")
  newqda <- matrixqda(c_mat, groups, priors, method = "t")
  newprior <- c(-1, 2)

  expect_error(
    predict(newlda, prior = newprior),
    "invalid 'prior'"
  )
  expect_error(
    predict(newqda, prior = newprior),
    "invalid 'prior'"
  )
  newprior <- c(.7, .7)
  expect_error(
    predict(newlda, prior = newprior),
    "invalid 'prior'"
  )
  expect_error(
    predict(newqda, prior = newprior),
    "invalid 'prior'"
  )

  expect_equal(sum(predict(newlda, newdata = matrix(
    0,
    nrow = 2, ncol = 2
  ))$posterior), 1)
  expect_equal(sum(predict(newlda, prior = c(.7, .3))$posterior[1, ]), 1)

  expect_equal(sum(predict(newqda, prior = c(.7, .3), newdata = matrix(
    0,
    nrow = 2, ncol = 2
  ))$posterior), 1)
  expect_equal(sum(predict(newqda)$posterior[1, ]), 1)
  expect_equal(class(logLik(ldamodel)), "logLik")
  expect_equal(class(logLik(qdamodel)), "logLik")
})

test_that("LDA/QDa_mat logLik works", {
  set.seed(20190628)
  ntotal <- 25
  covmatrix <- matrix(c(1, .5, .5, 1), nrow = 2)
  badcovmatrix <- matrix(c(1, .96, .96, 1), nrow = 2, ncol = 2)
  a_mat <- rmatrixnorm(ntotal,
    mean = matrix(0, nrow = 2, ncol = 2),
    U = covmatrix, V = covmatrix
  )
  b_mat <- rmatrixnorm(ntotal,
    mean = matrix(1, nrow = 2, ncol = 2),
    U = covmatrix, V = covmatrix
  )

  c_mat <- array(c(a_mat, b_mat), dim = c(2, 2, 2 * ntotal))
  groups <- c(rep(1, ntotal), rep(2, ntotal))
  priors <- c(.5, .5)

  # norm vs t DF 0
  ldamodel <- matrixlda(c_mat, groups, priors, method = "normal")
  qdamodel <- matrixqda(c_mat, groups, priors, method = "normal")
  ldamodelc <- matrixlda(c_mat, groups, priors, method = "t", nu = 0)
  qdamodelc <- matrixqda(c_mat, groups, priors, method = "t", nu = 0)

  expect_equal(ldamodel$method, ldamodelc$method)
  expect_equal(qdamodel$method, qdamodelc$method)

  # norm vs t
  ldamodel <- matrixlda(c_mat, groups, priors, method = "normal")
  qdamodel <- matrixqda(c_mat, groups, priors, method = "normal")
  ldamodelc <- matrixlda(c_mat, groups, priors, method = "t", nu = 5)
  qdamodelc <- matrixqda(c_mat, groups, priors, method = "t", nu = 5)

  expect_equal(
    attributes(logLik(ldamodel))$df,
    attributes(logLik(ldamodelc))$df
  )
  expect_equal(
    attributes(logLik(qdamodel))$df,
    attributes(logLik(qdamodelc))$df
  )

  # row.mean vs col.mean
  ldamodel <- matrixlda(c_mat, groups, priors,
    row.mean = TRUE,
    U = covmatrix, V = covmatrix
  )
  qdamodel <- matrixqda(c_mat, groups, priors,
    row.mean = TRUE,
    U = badcovmatrix, V = badcovmatrix
  )
  ldamodelc <- matrixlda(c_mat, groups, priors, col.mean = TRUE)
  qdamodelc <- matrixqda(c_mat, groups, priors, col.mean = TRUE)
  # only works because square!
  expect_equal(
    attributes(logLik(ldamodel))$df,
    attributes(logLik(ldamodelc))$df
  )
  expect_equal(
    attributes(logLik(qdamodel))$df,
    attributes(logLik(qdamodelc))$df
  )

  # row variance vs col variance AR
  ldamodel <- matrixlda(c_mat, groups, priors, row.variance = "AR")
  qdamodel <- matrixqda(c_mat, groups, priors, row.variance = "AR")
  ldamodelc <- matrixlda(c_mat, groups, priors, col.variance = "AR")
  qdamodelc <- matrixqda(c_mat, groups, priors, col.variance = "AR")
  # only works because square!
  expect_equal(
    attributes(logLik(ldamodel))$df,
    attributes(logLik(ldamodelc))$df
  )
  expect_equal(
    attributes(logLik(qdamodel))$df,
    attributes(logLik(qdamodelc))$df
  )
  # row variance vs col variance CS
  ldamodel <- matrixlda(c_mat, groups, priors, row.variance = "CS")
  qdamodel <- matrixqda(c_mat, groups, priors, row.variance = "CS")
  ldamodelc <- matrixlda(c_mat, groups, priors, col.variance = "CS")
  qdamodelc <- matrixqda(c_mat, groups, priors, col.variance = "CS")
  # only works because square!
  expect_equal(
    attributes(logLik(ldamodel))$df,
    attributes(logLik(ldamodelc))$df
  )
  expect_equal(
    attributes(logLik(qdamodel))$df,
    attributes(logLik(qdamodelc))$df
  )

  # row variance vs col variance corr
  ldamodel <- matrixlda(c_mat, groups, priors, row.variance = "corr")
  qdamodel <- matrixqda(c_mat, groups, priors, row.variance = "corr")
  ldamodelc <- matrixlda(c_mat, groups, priors, col.variance = "corr")
  qdamodelc <- matrixqda(c_mat, groups, priors, col.variance = "corr")
  # only works because square!
  expect_equal(
    attributes(logLik(ldamodel))$df,
    attributes(logLik(ldamodelc))$df
  )
  expect_equal(
    attributes(logLik(qdamodel))$df,
    attributes(logLik(qdamodelc))$df
  )


  # row variance vs col variance I
  ldamodel <- matrixlda(c_mat, groups, priors, row.variance = "I")
  qdamodel <- matrixqda(c_mat, groups, priors, row.variance = "I")
  ldamodelc <- matrixlda(c_mat, groups, priors, col.variance = "I")
  qdamodelc <- matrixqda(c_mat, groups, priors, col.variance = "I")
  # only works because square!
  expect_equal(
    attributes(logLik(ldamodel))$df,
    attributes(logLik(ldamodelc))$df
  )
  expect_equal(
    attributes(logLik(qdamodel))$df,
    attributes(logLik(qdamodelc))$df
  )

  # AR vs CS
  ldamodel <- matrixlda(c_mat, groups, priors, row.variance = "AR")
  qdamodel <- matrixqda(c_mat, groups, priors, row.variance = "AR")
  ldamodelc <- matrixlda(c_mat, groups, priors, row.variance = "CS")
  qdamodelc <- matrixqda(c_mat, groups, priors, row.variance = "CS")
  # only works because square!
  expect_equal(
    attributes(logLik(ldamodel))$df,
    attributes(logLik(ldamodelc))$df
  )
  expect_equal(
    attributes(logLik(qdamodel))$df,
    attributes(logLik(qdamodelc))$df
  )
  expect_lt(
    attributes(logLik(ldamodel))$df,
    attributes(logLik(qdamodel))$df
  )
})

test_that("Warning messages for inverted t", {
  expect_warning(dmatrixinvt(matrix(1, nrow = 5, ncol = 1),
    df = 5, V = 4
  ), "undefined")

  expect_warning(a_mat <- dmatrixinvt(t(matrix(1, nrow = 5, ncol = 1)),
    df = 5, U = 4
  ), "undefined")
  expect_true(is.nan(a_mat))
})


test_that("Sample size for ML:", {
  expect_error(MLmatrixt(rmatrixt(2, df = 5, mean = matrix(0, 5, 6))))

  expect_error(MLmatrixnorm(rmatrixnorm(2, mean = matrix(0, 5, 6))))
})
