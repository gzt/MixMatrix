library(matrixdist)

context("Testing input integrity")


test_that("trying wrong type of input", {
  expect_error(rmatrixnorm(n = 1, mean = matrix(c("A",1))))
  expect_error(rmatrixt(n = 1, df =1, mean = matrix(c("A",1))))
  expect_error(rmatrixinvt(n = 1, df =1, mean = matrix(c("A",1))))
  expect_error(dmatrixnorm(x = matrix(c("A",1))))
  expect_error(dmatrixt(x = matrix(c("A",1)), df =1))
  expect_error(dmatrixinvt(x = matrix(c("A",1)),df=1))

  expect_error(rmatrixnorm(n = 1, mean = matrix(c(0,0)), U = matrix(c("A",0,0,1),nrow=2)))
  expect_error(rmatrixt(n = 1, df = 1, mean = matrix(c(0,0)), U = matrix(c("A",0,0,1),nrow=2)))
  expect_error(rmatrixinvt(n = 1, df = 1, mean = matrix(c(0,0)), U = matrix(c("A",0,0,1),nrow=2)))
  expect_error(dmatrixnorm(x = matrix(c(0,0)), U = matrix(c("A",0,0,1),nrow=2)))
  expect_error(dmatrixt(x = matrix(c(0,0)), df = 1, U = matrix(c("A",0,0,1),nrow=2)))
  expect_error(dmatrixinvt(x =mean(matrix(c(0,0))), df =1, U = matrix(c("A",0,0,1),nrow=2)))


  expect_error(rmatrixnorm(n="A", mean = matrix(c(100,0,-100,0,25,-1000),nrow=2),
                  L= matrix(c(2,1,0,.1),nrow=2),
                  R = matrix(c(5,1,2,0,6,1,-1,2,10),nrow=3),list=FALSE))
  expect_error(rmatrixnorm(n = 10, mean = matrix(c(100,0,-100,0,25,-1000),nrow="A"),
                             L= matrix(c(2,1,0,.1),nrow=2),
                             R = matrix(c(5,1,2,0,6,1,-1,2,10),nrow=3),list=FALSE))
  expect_error(rmatrixnorm(n = 10, mean = matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             L= matrix(c("A",1,0,.1),nrow=2),
                             R = matrix(c(5,1,2,0,6,1,-1,2,10),nrow=3),list=FALSE))
  expect_error(rmatrixnorm(n = 10, mean = matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             L= matrix(c(2,1,0,.1),nrow=2),
                             R = matrix(c("A",1,2,0,6,1,-1,2,10),nrow=3),list=FALSE))
  expect_error(rmatrixnorm(n = 10, mean = matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             L= matrix(c(2,1,0,.1),nrow=2),
                             R = matrix(c("A",1,2,0,6,1,-1,2,10),nrow=3),list=TRUE))
  expect_error(rmatrixnorm(n = 10, mean = matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             U= matrix(c("A",1,0,.1),nrow=2),
                             R = matrix(c(5,1,2,0,6,1,-1,2,10),nrow=3),list=FALSE))

  expect_error(dmatrixnorm(matrix(c("A",0,-100,0,25,-1000),nrow=2),
              mean = matrix(c(100,0,-100,0,25,-1000),nrow=2),
              L= matrix(c(2,1,0,.1),nrow=2),log=TRUE ))
  expect_error(dmatrixnorm(matrix(c(100,0,-100,0,25,-1000),nrow="A"),
                             mean = matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             L= matrix(c(2,1,0,.1),nrow=2),log=TRUE ))
  expect_error(dmatrixnorm(matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             mean = matrix(c("A",0,-100,0,25,-1000),nrow=2),
                             L= matrix(c(2,1,0,.1),nrow=2),log=TRUE ))
  expect_error(dmatrixnorm(matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             mean = matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             L= matrix(c("A",1,0,.1),nrow=2),log=TRUE ))
  expect_error(dmatrixnorm(matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             mean = matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             L= matrix(c(2,1,0,.1),nrow="A"),log=TRUE ))
  expect_error(dmatrixnorm(matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             mean = matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             L= matrix(c(2,1,0,.1),nrow=2),
                             R = matrix(c("A",1,2,0,6,1,-1,2,10),nrow=3),log=TRUE ))
  expect_error(dmatrixnorm(matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             mean = matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             U= matrix(c("A",1,0,.1),nrow=2),log=TRUE ))
  expect_error(dmatrixnorm(matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             mean = matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             L= matrix(c(2,1,0,.1),nrow=2),
                             V = matrix(c("A",1,2,0,6,1,-1,2,10),nrow=3),log=TRUE ))


  expect_error(MLmatrixnorm(rmatrixnorm(n = 100, mean=diag(5)), tol="Q"))
  expect_error(MLmatrixnorm(rmatrixnorm(n = 100, mean=diag(5)), max.iter="Q"))
  expect_error(MLmatrixnorm(rmatrixnorm(n = 100, mean=diag(5)),
                              U = matrix("Q",nrow=5,ncol=5)))
  expect_error(MLmatrixnorm(rmatrixnorm(n = 100, mean=diag(5)),
                              V = matrix("Q",nrow=5,ncol=5)))

  expect_error(toepgenerate(7,"A"))
  expect_error(toepgenerate("A",.4))

  expect_error(rCholWishart(1,df=5,Sigma="A"))
  expect_error(rCholWishart(1,df="A",Sigma=diag(5)))
  expect_error(rCholWishart("A",df=5,Sigma=diag(5)))
  expect_error(rCholWishart(.4,df=5,Sigma=diag(5)))

  expect_error(rInvCholWishart(1,df=5,Sigma="A"))
  expect_error(rInvCholWishart(1,df="A",Sigma=diag(5)))
  expect_error(rInvCholWishart("A",df=4,Sigma=diag(5)))
  expect_error(rInvCholWishart(.4,df=4,Sigma=diag(5)))

})

test_that("Testing bad matrix dimension input",{
  A <- diag(3)
  B <- diag(4)

  expect_error(rmatrixnorm(n = 1, mean=A,L=B,R=A,list=FALSE))
  expect_error(rmatrixnorm(n = 1, mean=A,L=A,R=B,list=FALSE))
  expect_error(rmatrixnorm(n = 1, mean=A,U=B,R=A,list=FALSE))
  expect_error(rmatrixnorm(n = 1, mean=A,L=A,V=B,list=FALSE))
  expect_error(rmatrixnorm(n = 1, mean=A,U=A,V=B,list=FALSE))
  expect_error(rmatrixnorm(n = 1, mean=A,U=A,R=B,list=FALSE))

  expect_error(rmatrixt(n = 1, df =1, mean=A,L=B,R=A,list=FALSE))
  expect_error(rmatrixt(n = 1, df =1, mean=A,L=A,R=B,list=FALSE))
  expect_error(rmatrixt(n = 1, df =1, mean=A,U=B,R=A,list=FALSE))
  expect_error(rmatrixt(n = 1, df =1, mean=A,L=A,V=B,list=FALSE))
  expect_error(rmatrixt(n = 1, df =1, mean=A,U=A,V=B,list=FALSE))
  expect_error(rmatrixt(n = 1, df =1, mean=A,U=A,R=B,list=FALSE))

  expect_error(rmatrixinvt(n = 1, df =1, mean=A,L=B,R=A,list=FALSE))
  expect_error(rmatrixinvt(n = 1, df =1, mean=A,L=A,R=B,list=FALSE))
  expect_error(rmatrixinvt(n = 1, df =1, mean=A,U=B,R=A,list=FALSE))
  expect_error(rmatrixinvt(n = 1, df =1, mean=A,L=A,V=B,list=FALSE))
  expect_error(rmatrixinvt(n = 1, df =1, mean=A,U=A,V=B,list=FALSE))
  expect_error(rmatrixinvt(n = 1, df =1, mean=A,U=A,R=B,list=FALSE))

  expect_error(dmatrixnorm(x = A, mean=A,L=B,R=A))
  expect_error(dmatrixnorm(x = A, mean=A,L=A,R=B))
  expect_error(dmatrixnorm(x = A, mean=A,U=B,R=A))
  expect_error(dmatrixnorm(x = A, mean=A,L=A,V=B))
  expect_error(dmatrixnorm(x = A, mean=A,U=A,V=B))
  expect_error(dmatrixnorm(x = A, mean=A,U=A,R=B))

  expect_error(dmatrixt(x = A, df =1, mean=A,L=B,R=A))
  expect_error(dmatrixt(x = A, df =1, mean=A,L=A,R=B))
  expect_error(dmatrixt(x = A, df =1, mean=A,U=B,R=A))
  expect_error(dmatrixt(x = A, df =1, mean=A,L=A,V=B))
  expect_error(dmatrixt(x = A, df =1, mean=A,U=A,V=B))
  expect_error(dmatrixt(x = A, df =1, mean=A,U=A,R=B))

  expect_error(dmatrixinvt(x = A, df =1, mean=A,L=B,R=A))
  expect_error(dmatrixinvt(x = A, df =1, mean=A,L=A,R=B))
  expect_error(dmatrixinvt(x = A, df =1, mean=A,U=B,R=A))
  expect_error(dmatrixinvt(x = A, df =1, mean=A,L=A,V=B))
  expect_error(dmatrixinvt(x = A, df =1, mean=A,U=A,V=B))
  expect_error(dmatrixinvt(x = A, df =1, mean=A,U=A,R=B))

  expect_error(MLmatrixnorm(rmatrixnorm(n = 100, mean=A),U=B,V=A))
  expect_error(MLmatrixnorm(rmatrixnorm(n = 100, mean=A),U=A,V=B))

  expect_error(rCholWishart(1,10,matrix(c(3,1,1,1,1,3),nrow=2)))
  expect_error(rInvCholWishart(1,10,matrix(c(3,1,1,1,1,3),nrow=2)))


})

test_that("Out of bounds numeric input: ",{
  expect_error(lmvgamma(-1,5))
  expect_error(lmvgamma(1,-5))
  A <- diag(5)
  A[5,5] <- 0
  expect_error(rmatrixt(0, 1, matrix(0,nrow=5,ncol=5), U = diag(5), V = diag(5)))
  expect_error(rmatrixt(1, -1, matrix(0,nrow=5,ncol=5), U = diag(5), V = diag(5)))
  expect_error(rmatrixt(1, 1, matrix(0,nrow=5,ncol=5), U = -diag(5), V = diag(5)))
  expect_error(rmatrixt(1, 1, matrix(0,nrow=5,ncol=5), U = diag(5), V = -diag(5)))
  expect_error(rmatrixt(1, 1, matrix(0,nrow=5,ncol=5), L = A, V = diag(5)))
  expect_error(rmatrixt(1, 1, matrix(0,nrow=5,ncol=5), U = diag(5), R = A))
  expect_error(rmatrixt(1, 1, matrix(0,nrow=5,ncol=5), U = A, V = diag(5)))
  expect_error(rmatrixt(1, 1, matrix(0,nrow=5,ncol=5), U = diag(5), V = A))

  expect_error(rmatrixinvt(0, 1, matrix(0,nrow=5,ncol=5), U = diag(5), V = diag(5)))
  expect_error(rmatrixinvt(1, -1, matrix(0,nrow=5,ncol=5), U = diag(5), V = diag(5)))
  expect_error(rmatrixinvt(1, 1, matrix(0,nrow=5,ncol=5), U = -diag(5), V = diag(5)))
  expect_error(rmatrixinvt(1, 1, matrix(0,nrow=5,ncol=5), U = diag(5), V = -diag(5)))
  expect_error(rmatrixinvt(1, 1, matrix(0,nrow=5,ncol=5), L = A, V = diag(5)))
  expect_error(rmatrixinvt(1, 1, matrix(0,nrow=5,ncol=5), U = diag(5), R = A))
  expect_error(rmatrixinvt(1, 1, matrix(0,nrow=5,ncol=5), U = A, V = diag(5)))
  expect_error(rmatrixinvt(1, 1, matrix(0,nrow=5,ncol=5), U = diag(5), V = A))

  expect_error(rmatrixnorm(0, matrix(0,nrow=5,ncol=5), U = diag(5), V = diag(5)))
  expect_error(rmatrixnorm(1, matrix(0,nrow=5,ncol=5), U = -diag(5), V = diag(5)))
  expect_error(rmatrixnorm(1, matrix(0,nrow=5,ncol=5), U = diag(5), V = -diag(5)))
  expect_error(rmatrixnorm(1, matrix(0,nrow=5,ncol=5), L =A, V = diag(5)))
  expect_error(rmatrixnorm(1, matrix(0,nrow=5,ncol=5), U = diag(5), R = A))
  expect_error(rmatrixnorm(1, 1, matrix(0,nrow=5,ncol=5), L = A, V = diag(5)))
  expect_error(rmatrixnorm(1, 1, matrix(0,nrow=5,ncol=5), U = diag(5), V= A))

  expect_error(dmatrixt(1, matrix(0,nrow=5,ncol=5), U = diag(5), V = diag(5)))

  expect_error(dmatrixinvt(1, matrix(0,nrow=5,ncol=5), U = diag(5), V = diag(5)))

  expect_error(dmatrixt(df = -1, x = matrix(0,nrow=5,ncol=5), U = diag(5), V = diag(5)))
  expect_error(dmatrixt(df = 1, x = matrix(0,nrow=5,ncol=5), U = -diag(5), V = diag(5)))
  expect_error(dmatrixt(df = 1, x = matrix(0,nrow=5,ncol=5), U = diag(5), V = -diag(5)))
  expect_error(dmatrixt(df = 1, x = matrix(0,nrow=5,ncol=5), L = A, V = diag(5)))
  expect_error(dmatrixt(df = 1, x = matrix(0,nrow=5,ncol=5), U = diag(5), R = A))
  expect_error(dmatrixt(df = 1, x = matrix(0,nrow=5,ncol=5), U = A, V = diag(5)))
  expect_error(dmatrixt(df = 1, x = matrix(0,nrow=5,ncol=5), U = diag(5), V = A))

  expect_error(dmatrixinvt(df = -1, x = matrix(0,nrow=5,ncol=5), U = diag(5), V = diag(5)))
  expect_error(dmatrixinvt(df = 1, x = matrix(0,nrow=5,ncol=5), U = -diag(5), V = diag(5)))
  expect_error(dmatrixinvt(df = 1, x = matrix(0,nrow=5,ncol=5), U = diag(5), V = -diag(5)))
  expect_error(dmatrixinvt(df = 1, x = matrix(0,nrow=5,ncol=5), L = A, V = diag(5)))
  expect_error(dmatrixinvt(df = 1, x = matrix(0,nrow=5,ncol=5), U = diag(5), R = A))
  expect_error(dmatrixinvt(df = 1, x = matrix(0,nrow=5,ncol=5), U = A, V = diag(5)))
  expect_error(dmatrixinvt(df = 1, x = matrix(0,nrow=5,ncol=5), U = diag(5), V = A))

  expect_error(dmatrixnorm(matrix(0,nrow=5,ncol=5), U = -diag(5), V = diag(5)))
  expect_error(dmatrixnorm(matrix(0,nrow=5,ncol=5), U = diag(5), V = -diag(5)))
  expect_error(dmatrixnorm(matrix(0,nrow=5,ncol=5), L =A, V = diag(5)))
  expect_error(dmatrixnorm(matrix(0,nrow=5,ncol=5), U = diag(5), R = A))
  expect_error(dmatrixnorm(1, matrix(0,nrow=5,ncol=5), L = A, V = diag(5)))
  expect_error(dmatrixnorm(1, matrix(0,nrow=5,ncol=5), U = diag(5), V= A))


  expect_error(rCholWishart(1,10,-diag(5)))
  expect_error(rInvCholWishart(1,10,-diag(5)))


  expect_error(rCholWishart(1,4,diag(5)))
  expect_error(rInvCholWishart(1,4,diag(5)))


  expect_error(rCholWishart(-1,10,diag(5)))
  expect_error(rInvCholWishart(-1,10,diag(5)))

})

test_that("Bad rank in covariance:",{
  A = matrix(0,nrow=2,ncol=2)

  expect_error(rmatrixnorm(2, mean=matrix(c(0),nrow=2,ncol=2),
                           L = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_error(rmatrixnorm(2, mean=matrix(c(0),nrow=2,ncol=2),
                           R = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_length(rmatrixnorm(2, mean=matrix(c(0),nrow=2,ncol=2),
                             L = matrix(c(1,1,.5,.5),nrow=2,ncol=2), force = TRUE),8)
  expect_length(rmatrixnorm(2, mean=matrix(c(0),nrow=2,ncol=2),
                           R = matrix(c(1,1,.5,.5),nrow=2,ncol=2),force = TRUE),8)
  expect_error(rmatrixnorm(2, mean=matrix(c(0),nrow=2,ncol=2),
                           U = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_error(rmatrixnorm(2, mean=matrix(c(0),nrow=2,ncol=2),
                           V = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))

  expect_error(rmatrixt(2,5, mean=matrix(c(0),nrow=2,ncol=2),
                           L = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_error(rmatrixt(2, 5, mean=matrix(c(0),nrow=2,ncol=2),
                           R = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_length(rmatrixt(2, 5, mean=matrix(c(0),nrow=2,ncol=2),
                            R = matrix(c(1,1,.5,.5),nrow=2,ncol=2),force = TRUE),8)
  expect_error(rmatrixt(2, 5, mean=matrix(c(0),nrow=2,ncol=2),
                           U = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_error(rmatrixt(2, 5, mean=matrix(c(0),nrow=2,ncol=2),
                           V = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))


  expect_error(rmatrixinvt(2,5, mean=matrix(c(0),nrow=2,ncol=2),
                        L = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_error(rmatrixinvt(2, 5, mean=matrix(c(0),nrow=2,ncol=2),
                        R = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_error(rmatrixinvt(2, 5, mean=matrix(c(0),nrow=2,ncol=2),
                        U = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_error(rmatrixinvt(2, 5, mean=matrix(c(0),nrow=2,ncol=2),
                        V = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))

  expect_error(dmatrixnorm(A, mean=matrix(c(0),nrow=2,ncol=2),
                           L = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_error(dmatrixnorm(A, mean=matrix(c(0),nrow=2,ncol=2),
                           R = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_error(dmatrixnorm(A, mean=matrix(c(0),nrow=2,ncol=2),
                           U = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_error(dmatrixnorm(A, mean=matrix(c(0),nrow=2,ncol=2),
                           V = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))

  expect_error(dmatrixinvt(A,5,  mean=matrix(c(0),nrow=2,ncol=2),
                           L = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_error(dmatrixinvt(A, 5,  mean=matrix(c(0),nrow=2,ncol=2),
                           R = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_error(dmatrixinvt(A, 5,  mean=matrix(c(0),nrow=2,ncol=2),
                           U = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_error(dmatrixinvt(A, 5, mean=matrix(c(0),nrow=2,ncol=2),
                           V = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))

  expect_error(dmatrixt(A,5, mean=matrix(c(0),nrow=2,ncol=2),
                           L = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_error(dmatrixt(A, 5, mean=matrix(c(0),nrow=2,ncol=2),
                           R = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_error(dmatrixt(A, 5, mean=matrix(c(0),nrow=2,ncol=2),
                           U = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))
  expect_error(dmatrixt(A, 5, mean=matrix(c(0),nrow=2,ncol=2),
                           V = matrix(c(1,1,.5,.5),nrow=2,ncol=2)))


  expect_error(rCholWishart(1,10,matrix(c(1,1,1,1),nrow=2)))
  expect_error(rInvCholWishart(1,10,matrix(c(1,1,1,1),nrow=2)))



})
