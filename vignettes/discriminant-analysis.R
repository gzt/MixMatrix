## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "##"
)

## ----generatedata--------------------------------------------------------
set.seed(20180222)
library('matrixdist')
A <- rmatrixnorm(30, mean = matrix(0, nrow=2, ncol=3))
B <- rmatrixnorm(30, mean = matrix(c(1,0), nrow=2, ncol=3))
C <- rmatrixnorm(30, mean = matrix(c(0,1), nrow=2, ncol=3))
ABC <- array(c(A,B,C), dim = c(2,3,90))
groups <- factor(c(rep("A",30),rep("B",30),rep("C",30)))
prior = c(30,30,30)/90
matlda <- matrixlda(x = ABC,grouping = groups, prior = prior)
matqda <- matrixqda(x = ABC,grouping = groups, prior = prior)

## ----predict-------------------------------------------------------------
ABC[,,c(1,31,61)] # true class memberships: A, B, C

predict(matlda, ABC[,,c(1,31,61)])

predict(matqda,  ABC[,,c(1,31,61)])


## ----objectstructure-----------------------------------------------------
matlda

matqda


