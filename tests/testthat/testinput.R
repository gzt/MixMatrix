library(matrixdist)

context("Testing input integrity")


test_that("trying wrong type of input", {
  expect_error(rmatrixnorm(n=1,mean=matrix(c("A",1))))
  expect_error(rmatrixt(n=1,mean=matrix(c("A",1))))
  expect_error(rmatrixinvt(n=1,mean=matrix(c("A",1))))
  expect_error(dmatrixnorm(x=matrix(c("A",1))))
  expect_error(dmatrixt(x=matrix(c("A",1)),df=1))
  expect_error(dmatrixinvt(x=matrix(c("A",1))))

  expect_error(rmatrixnorm(n=1,mean=matrix(c(0,0)), U = matrix(c("A",0,0,1),nrow=2)))
  expect_error(rmatrixt(n=1,df = 1, mean=matrix(c(0,0)), U = matrix(c("A",0,0,1),nrow=2)))
  expect_error(rmatrixinvt(n=1,mean=matrix(c(0,0)), U = matrix(c("A",0,0,1),nrow=2)))
  expect_error(dmatrixnorm(x=matrix(c(0,0)), U = matrix(c("A",0,0,1),nrow=2)))
  expect_error(dmatrixt(x=matrix(c(0,0)),df = 1, U = matrix(c("A",0,0,1),nrow=2)))
  expect_error(dmatrixinvt(x=mean(matrix(c(0,0))),df=1, U = matrix(c("A",0,0,1),nrow=2)))


  expect_error(rmatrixnorm(n="A",mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
                  L=matrix(c(2,1,0,.1),nrow=2),
                  R = matrix(c(5,1,2,0,6,1,-1,2,10),nrow=3),list=FALSE))
  expect_error(rmatrixnorm(n=10,mean=matrix(c(100,0,-100,0,25,-1000),nrow="A"),
                             L=matrix(c(2,1,0,.1),nrow=2),
                             R = matrix(c(5,1,2,0,6,1,-1,2,10),nrow=3),list=FALSE))
  expect_error(rmatrixnorm(n=10,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             L=matrix(c("A",1,0,.1),nrow=2),
                             R = matrix(c(5,1,2,0,6,1,-1,2,10),nrow=3),list=FALSE))
  expect_error(rmatrixnorm(n=10,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             L=matrix(c(2,1,0,.1),nrow=2),
                             R = matrix(c("A",1,2,0,6,1,-1,2,10),nrow=3),list=FALSE))
  expect_error(rmatrixnorm(n=10,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             L=matrix(c(2,1,0,.1),nrow=2),
                             R = matrix(c("A",1,2,0,6,1,-1,2,10),nrow=3),list=TRUE))
  expect_error(rmatrixnorm(n=10,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             U=matrix(c("A",1,0,.1),nrow=2),
                             R = matrix(c(5,1,2,0,6,1,-1,2,10),nrow=3),list=FALSE))

  expect_error(dmatrixnorm(matrix(c("A",0,-100,0,25,-1000),nrow=2),
              mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
              L=matrix(c(2,1,0,.1),nrow=2),log=TRUE ))
  expect_error(dmatrixnorm(matrix(c(100,0,-100,0,25,-1000),nrow="A"),
                             mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             L=matrix(c(2,1,0,.1),nrow=2),log=TRUE ))
  expect_error(dmatrixnorm(matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             mean=matrix(c("A",0,-100,0,25,-1000),nrow=2),
                             L=matrix(c(2,1,0,.1),nrow=2),log=TRUE ))
  expect_error(dmatrixnorm(matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             L=matrix(c("A",1,0,.1),nrow=2),log=TRUE ))
  expect_error(dmatrixnorm(matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             L=matrix(c(2,1,0,.1),nrow="A"),log=TRUE ))
  expect_error(dmatrixnorm(matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             L=matrix(c(2,1,0,.1),nrow=2),
                             R = matrix(c("A",1,2,0,6,1,-1,2,10),nrow=3),log=TRUE ))
  expect_error(dmatrixnorm(matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             U=matrix(c("A",1,0,.1),nrow=2),log=TRUE ))
  expect_error(dmatrixnorm(matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
                             L=matrix(c(2,1,0,.1),nrow=2),
                             V = matrix(c("A",1,2,0,6,1,-1,2,10),nrow=3),log=TRUE ))


  expect_error(mle.matrixnorm(rmatrixnorm(n=100,mean=diag(5)), tol="Q"))
  expect_error(mle.matrixnorm(rmatrixnorm(n=100,mean=diag(5)), max.iter="Q"))
  expect_error(mle.matrixnorm(rmatrixnorm(n=100,mean=diag(5)),
                              U = matrix("Q",nrow=5,ncol=5)))
  expect_error(mle.matrixnorm(rmatrixnorm(n=100,mean=diag(5)),
                              V = matrix("Q",nrow=5,ncol=5)))

  expect_error(toepgenerate(7,"A"))
  expect_error(toepgenerate("A",.4))


})

test_that("Testing bad matrix dimension input",{
  A = diag(3)
  B = diag(4)

  expect_error(rmatrixnorm(n=1,mean=A,L=B,R=A,list=FALSE))
  expect_error(rmatrixnorm(n=1,mean=A,L=A,R=B,list=FALSE))
  expect_error(rmatrixnorm(n=1,mean=A,U=B,R=A,list=FALSE))
  expect_error(rmatrixnorm(n=1,mean=A,L=A,V=B,list=FALSE))
  expect_error(rmatrixnorm(n=1,mean=A,U=A,V=B,list=FALSE))
  expect_error(rmatrixnorm(n=1,mean=A,U=A,R=B,list=FALSE))

  expect_error(rmatrixt(n=1,df=1,mean=A,L=B,R=A,list=FALSE))
  expect_error(rmatrixt(n=1,df=1,mean=A,L=A,R=B,list=FALSE))
  expect_error(rmatrixt(n=1,df=1,mean=A,U=B,R=A,list=FALSE))
  expect_error(rmatrixt(n=1,df=1,mean=A,L=A,V=B,list=FALSE))
  expect_error(rmatrixt(n=1,df=1,mean=A,U=A,V=B,list=FALSE))
  expect_error(rmatrixt(n=1,df=1,mean=A,U=A,R=B,list=FALSE))

  expect_error(rmatrixinvt(n=1,df=1,mean=A,L=B,R=A,list=FALSE))
  expect_error(rmatrixinvt(n=1,df=1,mean=A,L=A,R=B,list=FALSE))
  expect_error(rmatrixinvt(n=1,df=1,mean=A,U=B,R=A,list=FALSE))
  expect_error(rmatrixinvt(n=1,df=1,mean=A,L=A,V=B,list=FALSE))
  expect_error(rmatrixinvt(n=1,df=1,mean=A,U=A,V=B,list=FALSE))
  expect_error(rmatrixinvt(n=1,df=1,mean=A,U=A,R=B,list=FALSE))

  expect_error(dmatrixnorm(x=A,mean=A,L=B,R=A))
  expect_error(dmatrixnorm(x=A,mean=A,L=A,R=B))
  expect_error(dmatrixnorm(x=A,mean=A,U=B,R=A))
  expect_error(dmatrixnorm(x=A,mean=A,L=A,V=B))
  expect_error(dmatrixnorm(x=A,mean=A,U=A,V=B))
  expect_error(dmatrixnorm(x=A,mean=A,U=A,R=B))

  expect_error(dmatrixt(x=A,df=1,mean=A,L=B,R=A))
  expect_error(dmatrixt(x=A,df=1,mean=A,L=A,R=B))
  expect_error(dmatrixt(x=A,df=1,mean=A,U=B,R=A))
  expect_error(dmatrixt(x=A,df=1,mean=A,L=A,V=B))
  expect_error(dmatrixt(x=A,df=1,mean=A,U=A,V=B))
  expect_error(dmatrixt(x=A,df=1,mean=A,U=A,R=B))

  expect_error(dmatrixinvt(x=A,df=1,mean=A,L=B,R=A))
  expect_error(dmatrixinvt(x=A,df=1,mean=A,L=A,R=B))
  expect_error(dmatrixinvt(x=A,df=1,mean=A,U=B,R=A))
  expect_error(dmatrixinvt(x=A,df=1,mean=A,L=A,V=B))
  expect_error(dmatrixinvt(x=A,df=1,mean=A,U=A,V=B))
  expect_error(dmatrixinvt(x=A,df=1,mean=A,U=A,R=B))

  expect_error(mle.matrixnorm(rmatrixnorm(n=100,mean=A),U=B,V=A))
  expect_error(mle.matrixnorm(rmatrixnorm(n=100,mean=A),U=A,V=B))

})
