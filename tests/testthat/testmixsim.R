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
