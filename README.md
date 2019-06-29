# MixMatrix

 [![Travis-CI Build Status](https://travis-ci.org/gzt/MixMatrix.svg?branch=master)](https://travis-ci.org/gzt/MixMatrix)
 [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/gzt/MixMatrix?branch=master&svg=true)](https://ci.appveyor.com/project/gzt/MixMatrix)
 [![Coverage Status](https://img.shields.io/codecov/c/github/gzt/MixMatrix/master.svg)](https://codecov.io/github/gzt/MixMatrix?branch=master)
[![packageversion](https://img.shields.io/badge/Package%20version-0.2.0%20-orange.svg?style=flat-square)](https://github.com/gzt/MixMatrix/releases)

A package for matrix variate distributions, including the normal, *t*, and 
inverted *t*. Currently LDA and QDA for matrix variate *t* distributions and 
normal distributions, EM for parameter estimation for matrix variate *t* 
distributions,  with sampling and density functions for those distributions 
as well as methods for parameter estimation for matrix variate normals and 
*t* with some restrictions on mean and variance parameters and some mixture 
modeling for matrix variate *t* distributions and normal distributions.

See the vignettes for an example of how it works.

It is currently possible to constrain the mean matrices for normal and *t* 
distributed matrices to have a common mean across rows, columns, or both, as 
well as AR(1), compound symmetric, on identity covariance matrices across rows, 
columns, or both in the maximum likelihood estimation functions and the LDA and
QDA functions, but not in the mixture models.


There isn't much R software out there for working with these distributions. The 
excellent [matrixsampling](https://cran.r-project.org/package=matrixsampling) 
package has sampling and distribution functions for these and many other matrix 
distributions,
as does [LaplacesDemon](https://cran.r-project.org/package=LaplacesDemon). The
[MatrixLDA](https://cran.r-project.org/package=MatrixLDA) package performs LDA 
for normal distributions with penalized likelihood.

The software can be installed by running:

    devtools::install_github("gzt/MixMatrix")
	
## Usage

The various `r*` and `d*` functions return or take as input an array indexed by 
the third dimension.

```
meanmatrix = matrix(1:12, nrow = 3)
A = rmatrixnorm(n = 10, mean = meanmatrix, U = diag(3), V = diag(4))
A[,,1:2]
dmatrixnorm(A, mean = meanmatrix, U = diag(3), V = diag(4), log = TRUE)
```

The package presents a method of maximum likelihood estimation of the parameters 
of the matrix variate *t* distribution using ECME.

```
X = rmatrixt(n = 100, mean = meanmatrix, U = diag(3), V = diag(4), df = 10)
MLmatrixt(X, fixed = FALSE) # fixed = FALSE indicates to estimate the DF parameter
```

Because it might be useful in conjunction with the `r*` and `d*` functions, the 
package also includes some convenience functions for generating AR(1) and 
compound symmetry correlation matrices.

```
ARgenerate(5, .5)
CSgenerate(5, .5)
```

The package also presents linear discriminant analysis and quadratic 
discriminant analysis for matrix variate distributions in the case of the 
normal and the *t*-distribution. In the case of the *t*, linear and quadratic 
refer whether the covariance matrices are constrained to be the same between 
classes rather than the form of the classifier.

```
 A <- rmatrixnorm(30,mean=matrix(0,nrow=3,ncol=4))
 B <- rmatrixnorm(30,mean=matrix(1,nrow=3,ncol=4), U = 2 * ARgenerate(3, .8))
 C <- array(c(A,B), dim=c(3,4,60))
 groups <- c(rep(1,30),rep(2,30))
 prior <- c(.5,.5)
 D <- matrixqda(C, groups, prior)
 predict(D)$posterior[1:10,]
```

In the future, this will include a more comprehensive treatment of matrix variate 
mixture modeling, including complete specification of covariance matrices in 
mixture modeling after the style of [Mclust](https://cran.r-project.org/package=mclust), 
[t-Eigen](https://cran.r-project.org/package=teigen), and other similar work.

## Contribution and contact information	

Please let me know if you have any issues or suggestions here: 
https://github.com/gzt/MixMatrix/issues

Please note that the `MixMatrix` project is released with a 
[Contributor Code of Conduct](https://gzt.github.io/MixMatrix/CODE_OF_CONDUCT.html). 
By contributing to this project, you agree to abide by its terms.

