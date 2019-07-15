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


    set.seed(20180221)
    A <- rmatrixnorm(30,mean=matrix(0,nrow=3,ncol=4))
    B <- rmatrixnorm(30,mean=matrix(1,nrow=3,ncol=4))
    C <- array(c(A,B), dim=c(3,4,60))
    prior <- c(.5,.5)

    rcmix <- (matrixmixture(C, prior = c(.5,.5), col.mean = TRUE, row.mean = TRUE))
    rmix <- (matrixmixture(C, prior = c(.5,.5), col.mean = FALSE, row.mean = TRUE))
    cmix <- (matrixmixture(C, prior = c(.5,.5), col.mean = TRUE, row.mean = FALSE))
    mix <- (matrixmixture(C, prior = c(.5,.5), col.mean = FALSE, row.mean = FALSE))

    
    trcmix <- (matrixmixture(C, prior = c(.5,.5), col.mean = TRUE, row.mean = TRUE, method = "t"))
    trmix <- (matrixmixture(C, prior = c(.5,.5), col.mean = FALSE, row.mean = TRUE, method = "t"))
    tcmix <- (matrixmixture(C, prior = c(.5,.5), col.mean = TRUE, row.mean = FALSE, method = "t"))
    tmix <- (matrixmixture(C, prior = c(.5,.5), col.mean = FALSE, row.mean = FALSE, method = "t"))
    

    llrcmix <- logLik(matrixmixture(C, prior = c(.5,.5), col.mean = TRUE, row.mean = TRUE))
    llrmix <- logLik(matrixmixture(C, prior = c(.5,.5), col.mean = FALSE, row.mean = TRUE))
    llcmix <- logLik(matrixmixture(C, prior = c(.5,.5), col.mean = TRUE, row.mean = FALSE))
    llmix <- logLik(matrixmixture(C, prior = c(.5,.5), col.mean = FALSE, row.mean = FALSE))

    
    lltrcmix <- logLik(matrixmixture(C, prior = c(.5,.5), col.mean = TRUE, row.mean = TRUE, method = "t"))
    lltrmix <- logLik(matrixmixture(C, prior = c(.5,.5), col.mean = FALSE, row.mean = TRUE, method = "t"))
    lltcmix <- logLik(matrixmixture(C, prior = c(.5,.5), col.mean = TRUE, row.mean = FALSE, method = "t"))
    lltmix <- logLik(matrixmixture(C, prior = c(.5,.5), col.mean = FALSE, row.mean = FALSE, method = "t"))



    expect_equal(rcmix$centers[1, 1, 1], rcmix$centers[1, 2, 1])
    expect_equal(rcmix$centers[1, 1, 1], rcmix$centers[2, 1, 1])

    expect_equal(trcmix$centers[1, 1, 1], trcmix$centers[1, 2, 1])
    expect_equal(trcmix$centers[1, 1, 1], trcmix$centers[2, 1, 1])
    

    expect_equal(rmix$centers[1, 1, 1], rmix$centers[1, 2, 1])
    expect_equal(cmix$centers[1, 1, 1], cmix$centers[2, 1, 1])

    expect_equal(trmix$centers[1, 1, 1], trmix$centers[1, 2, 1])
    expect_equal(tcmix$centers[1, 1, 1], tcmix$centers[2, 1, 1])
    
    expect_equal(attributes(llrcmix)$df, attributes(lltrcmix)$df)
    expect_equal(attributes(llmix)$df, attributes(lltmix)$df)
    expect_equal(attributes(llcmix)$df, attributes(lltcmix)$df)
    expect_equal(attributes(llrmix)$df, attributes(lltrmix)$df)
    expect_lt(attributes(llrcmix)$df,attributes(llcmix)$df)
    expect_lt(attributes(llcmix)$df,attributes(llmix)$df)
    expect_lt(attributes(llrmix)$df,attributes(llmix)$df)

    

})
