\name{simplex.sample}
\alias{simplex.sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample points uniformly distributed on a (constraint subset) of an n-dimensional simplex.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simplex.sample(N, n, constr = NULL, algo = "har")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{Number of points to generate}
  \item{n}{Dimensionality of the space}
  \item{constr}{Additional constraints (optional)}
  \item{algo}{Algorithm to be used (optional, default "har")}
}
\details{
The algorithm can be "har" (for Hit and Run), "bound" (for rejection sampling from the bounding box) or "simplex" (for rejection sampling from the simplex).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  \item{samples}{A matrix containing the generated samples}
  \item{miss}{The number of misses}
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
function (N, n, constr = NULL, algo = "har") 
{
    if (algo == "simplex") {
        hit <- function(x) {
            TRUE
        }
        if (!is.null(constr)) {
            hit <- function(x) {
                x <- c(x, -1)
                min(constr \%*\% x <= 0)
            }
        }
        result <- simplexReject(N, n, hit)
        samples <- t(sapply(result$samples, function(x) {
            x
        }))
        return(list(samples = samples, miss = result$miss))
    }
    else {
        basis <- simplex.basis(n)
        a <- simplex.createConstraints(basis, constr)
        hit <- simplex.createHit(a)
        bound <- createBoundBox(a)
        result <- NULL
        if (algo == "har") {
            result <- har(bound$start, N, bound$bound, hit)
        }
        else if (algo == "bound") {
            result <- boundingBoxReject(N, bound$lb, bound$ub, 
                hit)
        }
        samples <- simplex.transformResult(basis, result$samples)
        return(list(samples = samples, miss = result$miss))
    }
    stop(paste("Algorithm", algo, "not supported"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
