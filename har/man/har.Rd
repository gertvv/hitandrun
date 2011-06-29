\name{har}
\alias{har}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
"Hit and Run" sampler implementation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
har(x0, niter, bound, hit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x0}{Starting point (must be in the shape)}
  \item{niter}{Number of points to generate}
  \item{bound}{Bounding box function. Given a point x and a direction d, it will give two numbers, l <= 0 and u >= 0, that indicate that any point y > x + u * d or z < x + l * d is outside the target shape.}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x0, niter, bound, hit) 
{
    stopifnot(as.logical(hit(x0)))
    misses <- 0
    x <- as.list(rep(0, niter))
    n <- length(x0)
    x[[1]] <- x0
    for (i in 2:niter) {
        d <- randDir(n)
        bounds <- bound(x[[i - 1]], d)
        wasHit <- FALSE
        while (!wasHit) {
            l <- runif(1, bounds[1], bounds[2])[1]
            xN <- x[[i - 1]] + l * d
            if (hit(xN)) {
                x[[i]] <- xN
                wasHit <- TRUE
            }
            else {
                misses <- misses + 1
                if (l > 0) {
                  bounds[2] <- l
                }
                if (l < 0) {
                  bounds[1] <- l
                }
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
