# MixMatrix 0.2.7

* Fix Doxygen update to _PACKAGE issue
* Fix names in DESCRIPTION
* Fix "if (class(x) == ...)" issue.

# MixMatrix 0.2.6
* Update to fix M_PI issue with Rcpp update.
* Fix some non-HTTPS URLs.

# MixMatrix 0.2.5

* Update with DOI for associated paper.
* Fix an issue with AR(1) and CS covariances. The way these are specified 
means that in case of mis-specification, they may be non-monotonic if
they are mis-specified and results are near the the boundary. 
The functions were written such that it would not stop in this case. 

# MixMatrix 0.2.4

* Fixed a minor bug in DF calculation for mixture models.
* Fixed a minor bug in log likelihood calculation for mixture models.

# MixMatrix 0.2.2

* Added more citations and references.
* Added reference to the paper this package accompanies.
* A small series of bug fixes and internal function updates to more easily
  support restrictions on means and variances different functions.

# MixMatrix 0.2.1

* Added a number of citations and references to manuals and vignettes, including
  clarifying sources of any other functions.
* Added `logLik` functions for LDA, QDA, and MixMatrix objects

# MixMatrix 0.2.0

* Start to add mixture modeling - currently support unconstrained variances only
  and fixed `nu` parameters for the `model = "t"` option..
* Remove list support for matrix LDA and QDA, since they were just a pain
  to support for the `predict` methods. 


# MixMatrix 0.1.0

* Preparing for CRAN submission
* Add pkgdown (https://gzt.github.io/MixMatrix/)
* Add ORCID
* Minor documentation changes
* Fix `nu` bug for `matrixlda` and `matrixqda` when using normal distribution
* Add vignette for *t* distribution

# MixMatrix 0.0.0.9949

* changing name to MixMatrix

# matrixdist 0.0.0.9927

* include EM algorithm (actually an ECME) for estimation of parameters of 
  matrix-variate *t*-distributions
* include variance specifications for EM for t distribution
* include the possibility of restricting covariance matrices to a correlation 
  structure.

# matrixdist 0.0.0.9926

* split Wishart functions into CholWishart package
* migrate some internal functions to `RcppArmadillo` to improve speed
* add support for multiple observation input to `dmatrixnorm` function
* add support for multiple observation input to `dmatrixt` and `dmatrixinvt`
* dmatrix function now have large components in C++
* add error checking to `dmatrixinvt` - density only defined when matrix is 
  positive definite
* improve speed of `MLmatrixnorm` by migrating more to C++ and fixing some 
  R bottlenecks

# matrixdist 0.0.0.9925
 
* add multivariate `digamma` 
* add a couple tests for `matrixt`
* new plan: add fitting for t-distribution, LDA for t-distribution

# matrixdist 0.0.0.9924

* Wrote density functions for the Wishart and InvWishart distributions. 
  Broke the Wishart functions and `mvgamma` functions into a separate file 
  just for them as they are all related and I may break them into 
their own package.

# matrixdist 0.0.0.9923

* Rewrote `rInvCholWishart` to fix a bug where it was not actually a Cholesky 
  decomposition - while upper triangular,
`tcrossprod()` was InvWishart, not `crossprod()`. The method is a little slower 
	since it involves computing `rInvWishart` and then taking the Cholesky 
		decomposition, but it is faster than doing it in R.

# matrixdist 0.0.0.9922
 
* Added LDA and QDA functions with `predict()` methods. These functions
  are relatively fragile.
* Marked some functions as internal to clean up the NAMESPACE

# matrixdist 0.0.0.9921

* Added support for using IID structure for covariance matrices
* correction on restricted means - only true if known variance

# matrixdist 0.0.0.9920

* added change to mean estimation to account for restricted means
* added changes to documentation


# matrixdist 0.0.0.9919

* Added a `NEWS.md` file to track changes to the package.
* Changed the name of the maximum likelihood parameter fitting function to 
  `MLmatrixnorm` since `mle.---` could be confused for a method.



