library(MixMatrix)

context("Testing matrixmixture")

test_that("Testing bad input",{
    
    set.seed(20180221)
    A <- rmatrixnorm(30,mean=matrix(0,nrow=3,ncol=4))
    B <- rmatrixnorm(30,mean=matrix(2,nrow=3,ncol=4))
    C <- array(c(A,B), dim=c(3,4,60))
    prior <- c(.5,.5)
    init = list(centers = array(c(rep(0,12),rep(2,12)), dim = c(3,4,2)),
                U = array(c(diag(3), diag(3)), dim = c(3,3,2)),
                V = array(c(diag(4), diag(4)), dim = c(4,4,2))
                )
    expect_error(matrixmixture(C, init, prior = c(.1,.1)))
    expect_error(matrixmixture(C, init, prior = 0))
    expect_error(matrixmixture(C, init, prior = c(5,.1)))
    expect_error(matrixmixture(C, init, prior = c(-1,.1)))
    expect_error(matrixmixture(C, init))
    expect_error(matrixmixture(list(), prior = c(.5,.5),model="t",nu = 10))
       expect_error(matrixmixture(numeric(0), prior = c(.5,.5),model="t",nu = 10))

    
} )

test_that("Bad results warn or stop",{

    set.seed(20180221)
    A <- rmatrixnorm(30,mean=matrix(0,nrow=3,ncol=4))
    B <- rmatrixnorm(30,mean=matrix(2,nrow=3,ncol=4))
    C <- array(c(A,B), dim=c(3,4,60))
    prior <- c(.5,.5)
    init = list(centers = array(c(rep(0,12),rep(2,12)), dim = c(3,4,2)),
                U = array(c(diag(3), diag(3)), dim = c(3,3,2)),
                V = array(c(diag(4), diag(4)), dim = c(4,4,2))
                )
    expect_warning(capture.output(matrixmixture(C, init, prior = c(.5,.5), iter = 1, verbose = 100),
                                  type = "output"))
    expect_warning(matrixmixture(C, init, prior = 2,model="t",nu = 10, iter = 1))
    expect_warning(matrixmixture(C, K = 2,model="t",nu = 10, iter = 1))


    }
)

test_that("Mean restrictions work",{

    test_equal <- function(x, tolerance = 1e-6) all(abs(x - x[[1]]) < tolerance)

    set.seed(20180221)
    A <- rmatrixnorm(30,mean=matrix(0,nrow=3,ncol=4))
    B <- rmatrixnorm(30,mean=matrix(1,nrow=3,ncol=4))
    C <- array(c(A,B), dim=c(3,4,60))
    prior <- c(.5,.5)

    expect_true( test_equal(c(matrixmixture(C, prior = c(.5,.5), col.mean = TRUE, row.mean = TRUE)$centers[,,1])))
    expect_true( test_equal (c(matrixmixture(C, prior = c(.5,.5), col.mean = FALSE, row.mean = TRUE)$centers[1,,1])))
    expect_true( test_equal(matrixmixture(C, prior = c(.5,.5), col.mean = TRUE, row.mean = FALSE)$centers[,1,1]))
    expect_true( !test_equal(matrixmixture(C, prior = c(.5,.5), col.mean = FALSE, row.mean = FALSE)$centers[1,,1]))

    
    expect_true(test_equal(matrixmixture(C, prior = c(.5,.5), col.mean = TRUE, row.mean = TRUE, model = "t")$centers[,,1]))
    expect_true(test_equal(matrixmixture(C, prior = c(.5,.5), col.mean = FALSE, row.mean = TRUE, model = "t")$centers[1,,1]))
    expect_true(test_equal(matrixmixture(C, prior = c(.5,.5), col.mean = TRUE, row.mean = FALSE, model = "t")$centers[,1,1]))
    expect_true(!test_equal(matrixmixture(C, prior = c(.5,.5), col.mean = FALSE, row.mean = FALSE, model = "t")$centers[,1,1]))
    

    llrcmix <- logLik(matrixmixture(C, prior = c(.5,.5), col.mean = TRUE, row.mean = TRUE))
    llrmix <- logLik(matrixmixture(C, prior = c(.5,.5), col.mean = FALSE, row.mean = TRUE))
    llcmix <- logLik(matrixmixture(C, prior = c(.5,.5), col.mean = TRUE, row.mean = FALSE))
    llmix <- logLik(matrixmixture(C, prior = c(.5,.5), col.mean = FALSE, row.mean = FALSE))

    
    lltrcmix <- logLik(matrixmixture(C, prior = c(.5,.5), col.mean = TRUE, row.mean = TRUE, model = "t"))
    lltrmix <- logLik(matrixmixture(C, prior = c(.5,.5), col.mean = FALSE, row.mean = TRUE, model = "t"))
    lltcmix <- logLik(matrixmixture(C, prior = c(.5,.5), col.mean = TRUE, row.mean = FALSE, model = "t"))
    lltmix <- logLik(matrixmixture(C, prior = c(.5,.5), col.mean = FALSE, row.mean = FALSE, model = "t"))

    expect_equal(attributes(llrcmix)$df, attributes(lltrcmix)$df)
    expect_equal(attributes(llmix)$df, attributes(lltmix)$df)
    expect_equal(attributes(llcmix)$df, attributes(lltcmix)$df)
    expect_equal(attributes(llrmix)$df, attributes(lltrmix)$df)
    expect_lt(attributes(llrcmix)$df,attributes(llcmix)$df)
    expect_lt(attributes(llcmix)$df,attributes(llmix)$df)
    expect_lt(attributes(llrmix)$df,attributes(llmix)$df)

    

})


test_that("Predict Mix Model works", {


    set.seed(20180221)
    A <- rmatrixnorm(30,mean=matrix(0,nrow=3,ncol=4))
    B <- rmatrixnorm(30,mean=matrix(1,nrow=3,ncol=4))
    C <- array(c(A,B), dim=c(3,4,60))
    prior <- c(.5,.5)

    mix <- matrixmixture(C, prior = c(.5,.5))
    expect_error(predict(mix, newdata = matrix(0,nrow = 3, ncol = 2)),
                 "dimension")
    expect_error(predict(mix, newdata = (matrix(0,nrow = 2, ncol = 3))),
               "dimension")
    
    expect_equal(sum(predict(mix, newdata = matrix(
                                           0, nrow = 3, ncol = 4))$posterior), 1)
    expect_equal(sum(predict(mix, prior = c(.7,.3))$posterior[1,]), 1)

})

test_that("Init function works", {

    set.seed(20180221)
    A <- rmatrixnorm(30,mean=matrix(0,nrow=3,ncol=4))
    B <- rmatrixnorm(30,mean=matrix(1,nrow=3,ncol=4))
    C <- array(c(A,B), dim=c(3,4,60))
    prior <- c(.5,.5)

    testinit <- init_matrixmixture(C,K=2, centers = matrix(7,3,4), U = 4*diag(3), V = 3*diag(4))
    testinit_two <- init_matrixmixture(C,K=2, init=list(centers = matrix(7,3,4), U = 4*diag(3), V = 3*diag(4)))
    expect_equal(testinit$U[1,1,1],4)
    expect_equal(testinit$U[2,2,2],4)
    expect_equal(testinit$V[2,2,2],3)
    expect_equal(testinit$centers[1,1,2],7)
    expect_equal(testinit_two$U[1,1,1],4)
    expect_equal(testinit_two$U[2,2,2],4)
    expect_equal(testinit$V[2,2,2],3)
    expect_equal(testinit_two$centers[1,1,2],7)

    


})
