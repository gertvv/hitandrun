\name{har}
\alias{adaptiveHar}
\alias{nonAdaptiveHar}
\alias{gibbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
"Hit and Run" sampler implementation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
adaptiveHar(x0, niter, bound, hit)
nonAdaptiveHar(x0, niter, bound, hit)
gibbs(x0, niter, bound, hit)
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
