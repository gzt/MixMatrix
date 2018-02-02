\name{rmatrixnorm}
\alias{rmatrixnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rmatrixnorm(n, mean, L = diag(dim(mean)[1]), R = diag(dim(mean)[2]), U = L \%*\% t(L), V = t(R) \%*\% R, list = FALSE, array = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
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
  \item{list}{
%%     ~~Describe \code{list} here~~
}
  \item{array}{
%%     ~~Describe \code{array} here~~
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
function (n, mean, L = diag(dim(mean)[1]), R = diag(dim(mean)[2]), 
    U = L \%*\% t(L), V = t(R) \%*\% R, list = FALSE, array = NULL) 
{
    if (!(n > 0)) 
        stop("n must be > 0")
    dims = dim(mean)
    if (n == 1 && list == FALSE && is.null(array)) {
        return(rmatrixnorm.one(mean, U = U, V = V))
    }
    if (list) {
        return(lapply(1:n, function(x) rmatrixnorm.one(mean, 
            U = U, V = V)))
    }
    if (!(list) && !(is.null(array))) {
        if (!(array)) 
            warning("list FALSE and array FALSE, returning array")
    }
    result = array(data = NA, dim = c(n, dims), dimnames = list(NULL, 
        dimnames(mean)))
    for (i in 1:n) {
        result[i, , ] = rmatrixnorm.one(mean, U = U, V = V)
    }
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
