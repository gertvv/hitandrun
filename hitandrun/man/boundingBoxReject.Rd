\name{boundingBoxReject}
\alias{boundingBoxReject}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rejection sampling from a bounding box}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
boundingBoxReject(niter, lb, ub, hit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{niter}{Desired number of samples}
  \item{lb}{Vector of lower bounds (one for each dimension)}
  \item{ub}{Vector of upper bounds (one for each dimension)}
  \item{hit}{Hit calculation function. Given a point x, returns TRUE if it is within the target shape.}
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
\code{\link{constraint}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (niter, lb, ub, hit) 
{
    misses <- 0
    x <- as.list(rep(0, niter))
    n <- length(lb)
    d <- ub - lb
    for (i in 1:niter) {
        wasHit <- FALSE
        while (!wasHit) {
            xN <- lb + d * runif(n)
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
