library(MixMatrix)

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
  A <-    rmatrixnorm(      n = 1,      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
                L = matrix(c(2, 1, 0, .1), nrow = 2),
                list = FALSE
                )
  set.seed(2018020201)
  B <-    rmatrixnorm(      n = 10,      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
                L = matrix(c(2, 1, 0, .1), nrow = 2),      list = TRUE    )
  set.seed(2018020201)
  C <-      rmatrixnorm(      n = 10,      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
                  L = matrix(c(2, 1, 0, .1), nrow = 2),
                  list = FALSE)
  expect_equal(A, B[[1]])
  expect_equal(A, C[, , 1])
  expect_equal(dmatrixnorm(A,  mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2)),
               dmatrixnorm(B[[1]],  mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2)))
  expect_equal(dmatrixnorm(A,  mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2)),
               dmatrixnorm(C,  mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2))[1])
  
  set.seed(2018020202)
  A <-    rmatrixt(      n = 1,      df = 2,      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
                   L = matrix(c(2, 1, 0, .1), nrow = 2),
                   list = FALSE    )
  set.seed(2018020202)
  B <-    rmatrixt(      n = 1,      df = 2,      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
                   L = matrix(c(2, 1, 0, .1), nrow = 2),
                   list = TRUE    )
  set.seed(2018020202)
  C <-    rmatrixt(      n = 1,      df = 2,      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
      L = matrix(c(2, 1, 0, .1), nrow = 2),
      array = TRUE    )
  
  expect_equal(A, B[[1]])
  expect_equal(A, C[, , 1])
  expect_equal(dmatrixt(A, df = 2,  mean = matrix(c(100, 0, -100, 0, 25, -1000),  nrow = 2)),
               dmatrixt(B[[1]],df = 2,  mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2)))
    expect_equal(dmatrixt(A, df = 2,  mean = matrix(c(100, 0, -100, 0, 25, -1000),  nrow = 2)),
                 dmatrixt(C,df = 2,  mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2)))
  
  set.seed(2018020203)
  A <-    rmatrixinvt(      n = 1,      df = 2,      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
                      L = matrix(c(2, 1, 0, .1), nrow = 2),
                      list = FALSE    )
  set.seed(2018020203)
  B <-    rmatrixinvt(      n = 1,      df = 2,      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
      L = matrix(c(2, 1, 0, .1), nrow = 2),
      list = TRUE    )
  set.seed(2018020203)
  C <-    rmatrixinvt(      n = 1,      df = 2,
                      mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
                      L = matrix(c(2, 1, 0, .1), nrow = 2),
                      array = TRUE                      )
  expect_equal(A, B[[1]])
  expect_equal(A, C[, , 1])

  expect_equal( dmatrixnorm(array(1:10, dim = c(1,1,10)),0,U=4),
                dnorm(1:10,0,2))
  expect_equal( dmatrixnorm(array(1:10, dim = c(1,1,10)),0,V=16),
                dnorm(1:10,0,4))

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
  A = as.vector(rmatrixt(    1e4,    df = df,    0,    V = df,    list = FALSE  ))
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
  expect_equal(    dmatrixt(      x,      df,      U = U.one,      V = df * V.one,      log = T    ),
               dmatrixt(      t(x),      df,      U = U.two,      V = df * V.two,      log = T    ),
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
                U = CSgenerate(4,.2), V = ARgenerate(5, .8), list = TRUE)
  B <- MLmatrixnorm(A, row.mean = TRUE, V = ARgenerate(5,.8))
  expect_equal(B$U[1, 1], 1)
  expect_equal(B$V[1, 1], 1)

  expect_true(B$convergence)
  expect_equal(B$mean[1, 1], B$mean[1, 2])
 C <- MLmatrixnorm(A, col.mean = TRUE, row.mean = TRUE, U = CSgenerate(4,.2))
  expect_equal(C$mean[1, 1], C$mean[2, 1])
  expect_equal(C$mean[1, 1], C$mean[1, 2])
 

  C <- MLmatrixnorm(A, col.mean = TRUE)
  expect_equal(C$mean[1, 1], C$mean[2, 1])
  C <- MLmatrixnorm(A, row.variance = "CS")
  expect_equal(C$U[1,2],C$U[1,4])
  C <- MLmatrixnorm(A, col.variance = "CS")
  expect_equal(C$V[1,2],C$V[1,4])
  C <- MLmatrixnorm(A, row.variance = "AR(1)")
  expect_equal(C$U[2,1],C$U[3,2])
  C <- MLmatrixnorm(A, col.variance = "AR(1)")
  expect_equal(C$V[2,1],C$V[3,2])
  C <- MLmatrixnorm(A, row.variance = "corr")
  expect_equal(C$U[2,2],C$U[1,1])
  C <- MLmatrixnorm(A, col.variance = "corr")
  expect_equal(C$V[2,2],C$V[3,3])
  C <- MLmatrixnorm(A, row.variance = "I")
  expect_equal(C$U[1,2],0)
  C <- MLmatrixnorm(A, col.variance = "I")
  expect_equal(C$V[1,4],0)

  D <- MLmatrixt(A, col.mean = TRUE, U = CSgenerate(4,.2), V = ARgenerate(5,.8))
  expect_true(D$convergence)
  expect_warning(MLmatrixt(A, fixed = FALSE, max.iter = 2))
  expect_equal(D$U[1, 1], 1)
  expect_equal(D$V[1, 1], 1)
  expect_equal(D$mean[1, 1], D$mean[2, 1])
  C <- MLmatrixt(A, col.mean = TRUE, row.mean = TRUE)
  expect_equal(C$mean[1, 1], C$mean[2, 1])
  expect_equal(C$mean[1, 1], C$mean[1, 2])
   
  C <- MLmatrixt(A, col.variance = "CS")
  expect_equal(C$V[1,2],C$V[1,5])
  C <- MLmatrixt(A, row.variance = "CS")
  expect_equal(C$U[1,2],C$U[1,4])
  C <- MLmatrixt(A, row.variance = "AR(1)")
  expect_equal(C$U[2,1],C$U[3,2])
  C <- MLmatrixt(A, col.variance = "AR(1)")
  expect_equal(C$V[2,1],C$V[3,2])
  C <- MLmatrixt(A, col.variance = "corr")
  expect_equal(C$V[1,1],C$V[2,2])
  C <- MLmatrixt(A, row.variance = "corr")
  expect_equal(C$U[1,1],C$U[2,2])
  C <- MLmatrixt(A, row.variance = "I")
  expect_equal(C$U[1,2],0)
  C <- MLmatrixt(A, col.variance = "I")
  expect_equal(C$V[1,4],0)
})

context("Testing LDA/QDA output")
test_that("Output of LDA/QDA/Predict", {
    set.seed(20190628)
    A <- rmatrixnorm(4, mean = matrix(0, nrow = 2, ncol = 2))
    B <- rmatrixnorm(4, mean = matrix(1, nrow = 2, ncol = 2))
    set.seed(20190628)

    C <- array(c(A,B), dim = c(2,2,8))
    Czero <- C
    Czero[1,1,] <- 0
  D <- array(0, dim = c(2,2,4))
  E <- array(c(A,D), dim = c(2,2,8))
  groups <- c(rep(1,4),rep(2,4))
  groups.empty <- factor(rep("1",8), levels = c("1","2"))
  priors = c(.5,.5)
    ldamodel <- matrixlda(C, groups, priors, subset = rep(TRUE,8))
    qdamodel <- matrixqda(C, groups, priors, subset = rep(TRUE,8))
    expect_error(matrixlda(Czero, groups, priors, subset = rep(TRUE, 8)), "constant")
    expect_error(suppressWarnings(matrixqda(Czero, groups, priors, subset = rep(TRUE, 8))), "constant")
    expect_error(predict(ldamodel, newdata = matrix(0,nrow = 3, ncol = 2)),
                 "dimension")
    expect_error(predict(ldamodel, newdata = (matrix(0,nrow = 2, ncol = 3))),
               "dimension")
    expect_error(predict(qdamodel, newdata = matrix(0,nrow = 3, ncol = 2)),
                 "dimension")
    expect_error(predict(qdamodel, newdata = (matrix(0,nrow = 2, ncol = 3))),
                 "dimension")
    
    expect_equal(sum(predict(ldamodel, newdata = matrix(
                                           0, nrow = 2, ncol = 2))$posterior), 1)
    expect_equal(sum(predict(ldamodel, prior = c(.7,.3))$posterior[1,]), 1)
    
    expect_equal(sum(predict(qdamodel, prior=c(.7,.3),newdata = matrix(
                                                          0, nrow = 2, ncol = 2))$posterior), 1)
    expect_equal(sum(predict(qdamodel)$posterior[1,]), 1)

    llda = logLik(ldamodel)
    lqda = logLik(qdamodel)
    expect_equal(class(llda), "logLik")
    expect_equal(nobs(llda),8)
    expect_equal(class(lqda) , "logLik")
    expect_equal(nobs(lqda),8)
    expect_equal(attributes(llda)$nobs,nobs(ldamodel))
    expect_equal(attributes(lqda)$nobs,nobs(qdamodel))
  newlda <- matrixlda(C, groups, priors, method = "t")
  newqda <- matrixqda(C, groups, priors, method = "t")
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

    expect_equal(sum(predict(newlda, newdata = matrix(
                                         0, nrow = 2, ncol = 2))$posterior), 1)
    expect_equal(sum(predict(newlda, prior = c(.7,.3))$posterior[1,]), 1)
    
    expect_equal(sum(predict(newqda, prior=c(.7,.3),newdata = matrix(
                                                        0, nrow = 2, ncol = 2))$posterior), 1)
    expect_equal(sum(predict(newqda)$posterior[1,]), 1)
    expect_equal(class(logLik(ldamodel)) , "logLik")
    expect_equal(class(logLik(qdamodel)) , "logLik")

})

test_that("LDA/QDA logLik works",{
    set.seed(20190628)
    ntotal = 25
    covmatrix =  matrix(c(1,.5,.5,1),nrow=2)
    badcovmatrix = matrix(c(1,.96,.96,1), nrow=2, ncol=2)
  A <- rmatrixnorm(ntotal, mean = matrix(0, nrow = 2, ncol = 2), U =covmatrix, V = covmatrix)
  B <- rmatrixnorm(ntotal, mean = matrix(1, nrow = 2, ncol = 2), U= covmatrix, V = covmatrix)
  
  C <- array(c(A,B), dim = c(2,2,2*ntotal))
  groups <- c(rep(1,ntotal),rep(2,ntotal))
  priors = c(.5,.5)
    
  # norm vs t DF 0
  ldamodel <- matrixlda(C, groups, priors, method = "normal")
  qdamodel <- matrixqda(C, groups, priors, method = "normal")
  ldamodelc <- matrixlda(C, groups, priors, method = "t", nu  = 0)
  qdamodelc <- matrixqda(C, groups, priors, method = "t", nu = 0)
    
  expect_equal(ldamodel$method, ldamodelc$method)
  expect_equal(qdamodel$method, qdamodelc$method)
    
   # norm vs t 
  ldamodel <- matrixlda(C, groups, priors, method = "normal")
  qdamodel <- matrixqda(C, groups, priors, method = "normal")
  ldamodelc <- matrixlda(C, groups, priors, method = "t")
  qdamodelc <- matrixqda(C, groups, priors, method = "t")
    
  expect_equal(attributes(logLik(ldamodel))$df, attributes(logLik(ldamodelc))$df)
  expect_equal(attributes(logLik(qdamodel))$df, attributes(logLik(qdamodelc))$df) 

  # row.mean vs col.mean
  ldamodel <- matrixlda(C, groups, priors, row.mean = TRUE, U = covmatrix, V = covmatrix)
  qdamodel <- matrixqda(C, groups, priors, row.mean = TRUE, U = badcovmatrix, V = badcovmatrix)
  ldamodelc <- matrixlda(C, groups, priors, col.mean = TRUE)
  qdamodelc <- matrixqda(C, groups, priors, col.mean = TRUE)
  # only works because square!
  expect_equal(attributes(logLik(ldamodel))$df, attributes(logLik(ldamodelc))$df)
  expect_equal(attributes(logLik(qdamodel))$df, attributes(logLik(qdamodelc))$df)
  
  # row variance vs col variance AR
  ldamodel <- matrixlda(C, groups, priors, row.variance = "AR")
  qdamodel <- matrixqda(C, groups, priors, row.variance = "AR")
  ldamodelc <- matrixlda(C, groups, priors, col.variance = "AR")
  qdamodelc <- matrixqda(C, groups, priors, col.variance = "AR")
  # only works because square!
  expect_equal(attributes(logLik(ldamodel))$df, attributes(logLik(ldamodelc))$df)
  expect_equal(attributes(logLik(qdamodel))$df, attributes(logLik(qdamodelc))$df)
  # row variance vs col variance CS
  ldamodel <- matrixlda(C, groups, priors, row.variance = "CS")
  qdamodel <- matrixqda(C, groups, priors, row.variance = "CS")
  ldamodelc <- matrixlda(C, groups, priors, col.variance = "CS")
  qdamodelc <- matrixqda(C, groups, priors, col.variance = "CS")
  # only works because square!
  expect_equal(attributes(logLik(ldamodel))$df, attributes(logLik(ldamodelc))$df)
  expect_equal(attributes(logLik(qdamodel))$df, attributes(logLik(qdamodelc))$df)
  
  # row variance vs col variance corr
  ldamodel <- matrixlda(C, groups, priors, row.variance = "corr")
  qdamodel <- matrixqda(C, groups, priors, row.variance = "corr")
  ldamodelc <- matrixlda(C, groups, priors, col.variance = "corr")
  qdamodelc <- matrixqda(C, groups, priors, col.variance = "corr")
  # only works because square!
  expect_equal(attributes(logLik(ldamodel))$df, attributes(logLik(ldamodelc))$df)
  expect_equal(attributes(logLik(qdamodel))$df, attributes(logLik(qdamodelc))$df)

  
  # row variance vs col variance I
  ldamodel <- matrixlda(C, groups, priors, row.variance = "I")
  qdamodel <- matrixqda(C, groups, priors, row.variance = "I")
  ldamodelc <- matrixlda(C, groups, priors, col.variance = "I")
  qdamodelc <- matrixqda(C, groups, priors, col.variance = "I")
  # only works because square!
  expect_equal(attributes(logLik(ldamodel))$df, attributes(logLik(ldamodelc))$df)
  expect_equal(attributes(logLik(qdamodel))$df, attributes(logLik(qdamodelc))$df)

  # AR vs CS
  ldamodel <- matrixlda(C, groups, priors, row.variance = "AR")
  qdamodel <- matrixqda(C, groups, priors, row.variance = "AR")
  ldamodelc <- matrixlda(C, groups, priors, row.variance = "CS")
  qdamodelc <- matrixqda(C, groups, priors, row.variance = "CS")
  # only works because square!
  expect_equal(attributes(logLik(ldamodel))$df, attributes(logLik(ldamodelc))$df)
  expect_equal(attributes(logLik(qdamodel))$df, attributes(logLik(qdamodelc))$df)
  expect_lt(attributes(logLik(ldamodel))$df, attributes(logLik(qdamodel))$df)

})

test_that("Warning messages for inverted t",{

  expect_warning(dmatrixinvt(matrix(1, nrow=5, ncol=1), df = 5, V = 4), "undefined")

  expect_warning(A <- dmatrixinvt(t(matrix(1, nrow=5, ncol=1)), df = 5, U = 4), "undefined")
  # A <- dmatrixinvt(matrix(1, nrow=5, ncol=1), df = 5, V = 4)
  expect_true(is.nan(A))



})


test_that("Sample size for ML:",{
    
    expect_error(MLmatrixt(rmatrixt(2,df = 5, mean = matrix(0,5,6))))
    
    expect_error(MLmatrixnorm(rmatrixnorm(2,mean = matrix(0,5,6))))

})
