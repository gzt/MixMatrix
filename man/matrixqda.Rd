% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matnormlda.R
\name{matrixqda}
\alias{matrixqda}
\title{Quadratic Discriminant Analysis for Matrix Variate Observations}
\usage{
matrixqda(x, grouping, prior, tol = 1e-04, ...)
}
\arguments{
\item{x}{3-D array or list of matrix data.}

\item{grouping}{vector}

\item{prior}{a vector of prior probabilities of the same length
as the number of classes}

\item{tol}{by default, \code{1e-4}. Tolerance parameter checks
for 0 variance.}

\item{...}{Arguments passed to or from other methods, such
as additional parameters to pass to \code{MLmatrixnorm} (e.g.,
\code{row.mean})}
}
\value{
Returns a list of class \code{matrixqda} containing
   the following components:
   \describe{
      \item{\code{prior}}{the prior probabilities used.}
      \item{\code{counts}}{the counts of group membership}
      \item{\code{means}}{the group means.}
      \item{\code{U}}{the between-row covariance matrices}
      \item{\code{V}}{the between-column covariance matrices}
      \item{\code{lev}}{levels of the grouping factor}
      \item{\code{N}}{The number of observations used.}
      \item{\code{call}}{The (matched) function call.}
   }
}
\description{
See \code{matrixlda}: quadratic discriminant analysis for matrix
variate observations.
}
\details{
This uses \code{MLmatrixnorm} to find the means and variances.
}
\examples{
#' set.seed(20180221)
A <- rmatrixnorm(30,mean=matrix(0,nrow=3,ncol=4))
B <- rmatrixnorm(30,mean=matrix(1,nrow=3,ncol=4))
C <- array(c(A,B), dim=c(3,4,60))
groups <- c(rep(1,30),rep(2,30))
prior <- c(.5,.5)
D <- matrixqda(C, groups, prior)
}
\seealso{
Other matrixlda: \code{\link{matrixlda}},
  \code{\link{predict.matrixlda}},
  \code{\link{predict.matrixqda}}
}