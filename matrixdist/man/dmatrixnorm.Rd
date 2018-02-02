\name{dmatrixnorm}
\alias{dmatrixnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dmatrixnorm(x, mean = array(0L, dim(x)), L = diag(dim(x)[1]), R = diag(dim(x)[2]), U = L \%*\% t(L), V = t(R) \%*\% R, log = FALSE)
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
    p = dim(U)[1]
    n = dim(V)[1]
    detU = det(U)
    detV = det(V)
    if (!(detU > 0 || detV > 0)) 
        stop("non-invertible matrix")
    Uinv = solve(U)
    Vinv = solve(V)
    XM = x - mean
    logresult = -0.5 * n * p * log(2 * pi) - 0.5 * n * log(detU) - 
        0.5 * p * log(detV) - 0.5 * sum(diag(Vinv \%*\% t(XM) \%*\% 
        Uinv \%*\% (XM)))
    if (log) 
        return(logresult)
    else return(exp(logresult))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
