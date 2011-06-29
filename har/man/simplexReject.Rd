\name{simplexReject}
\alias{simplexReject}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs rejection sampling within an n-dimensional simplex.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simplexReject(niter, n, hit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{niter}{Number of points to generate}
  \item{n}{Dimensionality of the space}
  \item{hit}{Hit calculation function. Given a point x, returns TRUE if it is within the target shape. Points are guaranteed to be within the simplex, so that does not need to be checked.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{samples}{A matrix containing the generated samples}
  \item{miss}{The number of misses}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Gert van Valkenhoef
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
function (niter, n, hit) 
{
    misses <- 0
    x <- as.list(rep(0, niter))
    for (i in 1:niter) {
        wasHit <- FALSE
        while (!wasHit) {
            xN <- simplex.baseSample(n, 1)
            if (hit(xN)) {
                x[[i]] <- xN
                wasHit <- TRUE
            }
            else {
                misses <- misses + 1
            }
        }
    }
    list(samples = x, miss = misses)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
