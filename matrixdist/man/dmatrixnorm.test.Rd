\name{dmatrixnorm.test}
\alias{dmatrixnorm.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dmatrixnorm.test(x, mean = array(0L, dim(x)), L = diag(dim(x)[1]), R = diag(dim(x)[2]), U = L \%*\% t(L), V = t(R) \%*\% R, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{mean}{
%%     ~~Describe \code{mean} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{R}{
%%     ~~Describe \code{R} here~~
}
  \item{U}{
%%     ~~Describe \code{U} here~~
}
  \item{V}{
%%     ~~Describe \code{V} here~~
}
  \item{log}{
%%     ~~Describe \code{log} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, mean = array(0L, dim(x)), L = diag(dim(x)[1]), R = diag(dim(x)[2]), 
    U = L \%*\% t(L), V = t(R) \%*\% R, log = FALSE) 
{
    dims = dim(x)
    if (!(dims[1] == dim(U)[2] && dim(U)[1] == dim(U)[2] && dims[2] == 
        dim(V)[1] && dim(V)[1] == dim(V)[2])) {
        stop("Non-conforming dimensions.")
    }
    vecx = as.vector(x)
    meanx = as.vector(mean)
    VU = V \%x\% U
    p = dim(U)[1]
    n = dim(V)[1]
    detVU = det(VU)
    if (!(detVU > 0)) 
        stop("non-invertible matrix")
    UVinv = solve(VU)
    XM = vecx - meanx
    logresult = -0.5 * n * p * log(2 * pi) - 0.5 * log(det(VU)) - 
        0.5 * sum(diag(t(XM) \%*\% UVinv \%*\% XM))
    if (log) 
        return(logresult)
    else return(exp(logresult))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
