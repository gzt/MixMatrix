## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,comment="")
library(matrixdist)

## ----mnorm---------------------------------------------------------------
set.seed(20180203)
x <- matrix(rnorm(6),nrow=3)
x
set.seed(20180203)
y <- rmatrixnorm(n = 1, mean=matrix(rep(0,6),nrow=3))
y
U <- 5 * diag(3) + 1
V <- matrix(c(2,0,0,.1),nrow=2)
mu = matrix(1:6,nrow=3)
set.seed(20180203)
z <- rmatrixnorm(n = 1, mean=mu,U=U,V=V)
z
mu + t(chol(U)) %*% y %*% chol(V)


## ----noptions------------------------------------------------------------
set.seed(20180203)
x <- rmatrixnorm(n = 1, mean=matrix(rep(0,6),nrow=3))
x
set.seed(20180203)
y <- rmatrixnorm(n = 100, mean=matrix(rep(0,6),nrow=3),list = TRUE)
y[[1]]
set.seed(20180203)
z <- rmatrixnorm(n = 100, mean=matrix(rep(0,6),nrow=3),array = TRUE)
z[ , , 1]


## ----density-------------------------------------------------------------
set.seed(20180202)
A = rmatrixnorm(n=1,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
L=matrix(c(2,1,0,.1),nrow=2))
dmatrixnorm(A,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
  L=matrix(c(2,1,0,.1),nrow=2),log=TRUE )


## ----mleone, cache=T-----------------------------------------------------
set.seed(20180202)
A = rmatrixnorm(n=100,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
   L=matrix(c(2,1,0,.1),nrow=2),list=TRUE)
results=MLmatrixnorm(A)
print(results)

## ----mlerow, cache = T---------------------------------------------------
results.fixrows = MLmatrixnorm(A, row.mean = TRUE, max.iter = 5)
print(results.fixrows)
# this failure is expected with misspecification! The number of iterations is also 
# fixed to be low so the vignette compiles quickly.

## ----mlevar,warning=FALSE------------------------------------------------
# tolerance and maximum iterations are limited here to make the vignette compile faster

B = rmatrixnorm(n = 50, mean = matrix(1:15,nrow = 5), 
                U = 3*stats::toeplitz(.6^(1:5)))
MLmatrixnorm(B, row.variance = "AR(1)", tol = 1e-5)
MLmatrixnorm(B, tol = 1e-5)


