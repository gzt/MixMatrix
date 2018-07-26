# matrixdist

 [![Travis-CI Build Status](https://travis-ci.org/gzt/matrixdist.svg?branch=master)](https://travis-ci.org/gzt/matrixdist)
 [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/gzt/matrixdist?branch=master&svg=true)](https://ci.appveyor.com/project/gzt/matrixdist)
 [![Coverage Status](https://img.shields.io/codecov/c/github/gzt/matrixdist/master.svg)](https://codecov.io/github/gzt/matrixdist?branch=master)

A package for matrix variate distributions, including the normal, *t*, and inverted *t*. 
Currently with sampling and density functions for those distributions as well as methods for 
parameter estimation for matrix variate normals with some restrictions on mean and variance
parameters. This also has EM for parameter estimation for matrix variate *t*-distributions, which is new.

See the [vignette](../vignettes/matrixnormal.html) for more.

I plan to expand this to include introduce some further possible restrictions on variance matrices. The $t$-distributions do not, at the moment, have any such restrictions implemented.


There isn't much R software out there for working with these distributions. The 
excellent [matrixsampling](https://cran.r-project.org/package=matrixsampling) package 
has sampling and distribution functions for these and many other matrix distributions,
as does [LaplacesDemon](https://cran.r-project.org/package=LaplacesDemon).

The software can be installed by running:

    devtools::install_github("gzt/matrixdist")

Please let me know if you have any issues or suggestions here: https://github.com/gzt/matrixdist/issues


