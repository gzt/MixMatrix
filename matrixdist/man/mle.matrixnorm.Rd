\name{mle.matrixnorm}
\alias{mle.matrixnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mle.matrixnorm(data, row.restrict = "none", col.restrict = "none", tol = 1e-09, max.iter = 100, U, V)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{row.restrict}{
%%     ~~Describe \code{row.restrict} here~~
}
  \item{col.restrict}{
%%     ~~Describe \code{col.restrict} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{max.iter}{
%%     ~~Describe \code{max.iter} here~~
}
  \item{U}{
%%     ~~Describe \code{U} here~~
}
  \item{V}{
%%     ~~Describe \code{V} here~~
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
function (data, row.restrict = "none", col.restrict = "none", 
    tol = 1e-09, max.iter = 100, U, V) 
{
    if (class(data) == "list") 
        data = aperm(array(unlist(data), dim = c(nrow(data[[1]]), 
            ncol(data[[1]]), length(data))), perm = c(3, 1, 2))
    dims = dim(data)
    if (max(dims[2]/dims[3], dims[3]/dims[2]) < (dims[1] - 1)) 
        warning("Need more observations to estimate parameters.")
    if (missing(U)) 
        U = diag(dims[2])
    if (missing(V)) 
        V = diag(dims[3])
    mu = apply(data, c(2, 3), mean)
    swept.data = sweep(data, c(2, 3), mu)
    iter = 0
    error.term = 1e+40
    while (iter < max.iter && error.term > tol) {
        inter.V = apply(swept.data, 1, function(x) (t(x) \%*\% 
            solve(U) \%*\% x))
        new.V = matrix(apply(inter.V, 1, sum), nrow = dims[3])/(dims[1] * 
            dims[2])
        inter.U = apply(swept.data, 1, function(x) ((x) \%*\% solve(V) \%*\% 
            t(x)))
        new.U = matrix(apply(inter.U, 1, sum), nrow = dims[2])/(dims[1] * 
            dims[3])
        new.U = new.U/(new.U[1, 1])
        error.term = max(sum((new.V - V)^2), sum((new.U - U)^2))
        V = new.V
        U = new.U
        iter = iter + 1
    }
    if (iter >= max.iter || error.term > tol) 
        warning("Failed to converge")
    return(list(mean = mu, U = U, V = V, iter = iter, tol = error.term))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
