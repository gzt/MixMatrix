% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mlmatrixt.R
\name{MLmatrixt}
\alias{MLmatrixt}
\title{Maximum likelihood estimation for matrix variate t distributions}
\usage{
MLmatrixt(
  data,
  row.mean = FALSE,
  col.mean = FALSE,
  row.variance = "none",
  col.variance = "none",
  df = 10,
  fixed = TRUE,
  tol = .Machine$double.eps^0.5,
  max.iter = 5000,
  U,
  V,
  ...
)
}
\arguments{
\item{data}{Either a list of matrices or a 3-D array with matrices in
dimensions 1 and 2, indexed by dimension 3.}

\item{row.mean}{By default, \code{FALSE}. If \code{TRUE}, will fit a
common mean within each row. If both this and \code{col.mean} are
\code{TRUE}, there will be a common mean for the entire matrix.}

\item{col.mean}{By default, \code{FALSE}. If \code{TRUE}, will fit a
common mean within each row. If both this and \code{row.mean} are
\code{TRUE}, there will be a common mean for the entire matrix.}

\item{row.variance}{Imposes a variance structure on the rows. Either
'none', 'AR(1)', 'CS' for 'compound symmetry', 'Correlation' for a
correlation matrix, or 'Independence' for
independent and identical variance across the rows.
Only positive correlations are allowed for AR(1) and CS and these
restrictions may not be guaranteed to converge.
Note that while maximum likelihood estimators are available (and used)
for the unconstrained variance matrices, \code{optim} is used for any
constraints so it may be considerably slower.}

\item{col.variance}{Imposes a variance structure on the columns.
Either 'none', 'AR(1)', 'CS', 'Correlation', or 'Independence'.
Only positive correlations are allowed for
AR(1) and CS.}

\item{df}{Starting value for the degrees of freedom. If \code{fixed = TRUE},
then this is required and not updated. By default, set to 10.}

\item{fixed}{Whether \code{df} is estimated or fixed.
By default, \code{TRUE}.}

\item{tol}{Convergence criterion. Measured against square deviation
between iterations of the two variance-covariance matrices.}

\item{max.iter}{Maximum possible iterations of the algorithm.}

\item{U}{(optional) Can provide a starting point for the U matrix.
By default, an identity matrix.}

\item{V}{(optional) Can provide a starting point for the V matrix.
By default, an identity matrix.}

\item{...}{(optional) additional arguments can be passed to \code{optim}
if using restrictions on the variance.}
}
\value{
Returns a list with the following elements:
\describe{
\item{\code{mean}}{the mean matrix}
\item{\code{U}}{the between-row covariance matrix}
\item{\code{V}}{the between-column covariance matrix}
\item{\code{var}}{the scalar variance parameter
(the first entry of the covariances are restricted to unity)}
\item{\code{nu}}{the degrees of freedom parameter}
\item{\code{iter}}{the number of iterations}
\item{\code{tol}}{the squared difference between iterations of
the variance matrices at the time of stopping}
\item{\code{logLik}}{log likelihood of result.}
\item{\code{convergence}}{a convergence flag,
\code{TRUE} if converged.}
\item{\code{call}}{The (matched) function call.}
}
}
\description{
For the matrix variate normal distribution, maximum likelihood estimates
exist for \eqn{N > max(p/q,q/p)+1} and are unique for \eqn{N > max(p,q)}.
The number necessary for the matrix variate t has not been worked out but
this is a lower bound. This implements an ECME algorithm to estimate the
mean, covariance, and degrees of freedom parameters. An AR(1), compound
symmetry, or independence restriction can be proposed for either or both
variance matrices. However, if they are inappropriate for the data, they may
fail with a warning.
}
\examples{
set.seed(20180202)
# drawing from a distribution with specified mean and covariance
A <- rmatrixt(
  n = 100, mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
  L = matrix(c(2, 1, 0, .1), nrow = 2), list = TRUE, df = 5
)
# fitting maximum likelihood estimates
results <- MLmatrixt(A, tol = 1e-5, df = 5)
print(results)
}
\references{
\if{html}{\out{<div class="sourceCode">}}\preformatted{Thompson, G Z.  R Maitra, W Q Meeker, A Bastawros (2019),
"Classification with the matrix-variate-t distribution", arXiv
e-prints arXiv:1907.09565 <https://arxiv.org/abs/1907.09565>

Dickey, James M. 1967. “Matricvariate Generalizations of the
Multivariate t Distribution and the Inverted Multivariate t
Distribution.” Ann. Math. Statist. 38 (2): 511–18.
\doi{10.1214/aoms/1177698967}

Liu, Chuanhai, and Donald B. Rubin. 1994. “The ECME Algorithm:
A Simple Extension of EM and ECM with Faster Monotone Convergence.”
Biometrika 81 (4): 633–48.
      \doi{10.2307/2337067}
}\if{html}{\out{</div>}}

Meng, Xiao-Li, and Donald B. Rubin. 1993. “Maximum Likelihood Estimation
via the ECM Algorithm: A General Framework.” Biometrika 80 (2): 267–78.
\doi{10.1093/biomet/80.2.267}

\if{html}{\out{<div class="sourceCode">}}\preformatted{Rubin, D.B. 1983. “Encyclopedia of Statistical Sciences.” In, 4th ed.,
  272–5. John Wiley.
}\if{html}{\out{</div>}}
}
\seealso{
\code{\link[=rmatrixnorm]{rmatrixnorm()}}, \code{\link[=rmatrixt]{rmatrixt()}},
\code{\link[=MLmatrixnorm]{MLmatrixnorm()}}
}
