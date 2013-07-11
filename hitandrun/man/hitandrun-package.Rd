\name{hitandrun-package}
\alias{hitandrun-package}
\alias{hitandrun}
\docType{package}
\title{
"Hit and Run" sampling
}
\description{
This package provides a "Hit and Run" sampler that generates a Markov chain whose stable state converges on the uniform distribution over a convex polytope.
The polytope is given by a set of constraints in standard linear programming form (Ax <= b).

Utilities are provided for sampling from subsets of the unit simplex (i.e. random variates that can be interpreted as weights satisfying certain constraints) and for specifying common constraints.
}
\details{
\code{\link{har}} is the core "Hit and Run" sampler, \code{\link{bbReject}} is the bounding box rejection sampler and \code{\link{simplex.sample}} samples uniformly from the unit simplex.

See \code{\link{simplex.createTransform}} and \code{\link{simplex.createConstraints}} for sampling from subsets of the unit simplex.
Utilities to specify common constraints are described in \code{\link{harConstraints}}.
}
\author{
Maintainer: Gert van Valkenhoef <g.h.m.van.valkenhoef@rug.nl>
}
\references{
    Tervonen, T., van Valkenhoef, G., Basturk, N., and Postmus, D. (2012) "Hit-And-Run enables efficient weight generation for simulation-based multiple criteria decision analysis". \emph{European Journal of Operational Research} 224(3) 552-559. \href{http://dx.doi.org/10.1016/j.ejor.2012.08.026}{doi:10.1016/j.ejor.2012.08.026}
}
\note{
"Hit and Run" is a Markov Chain Monte Carlo (MCMC) method, so generated samples form a correlated time series. To get a uniform sample, you need O*(n^3) samples, where n is the dimension of the sampling space.
}
\keyword{hit-and-run}
\seealso{
\code{\link{har}}
\code{\link{bbReject}}
\code{\link{simplex.sample}}

\code{\link{simplex.createTransform}}
\code{\link{simplex.createConstraints}}

\code{\link{harConstraints}}

\code{\link{createSeedPoint}}
\code{\link{createBoundBox}}
}
\examples{
# Example: sample weight vectors where w_1 >= w_2 and w_1 >= w_3
n <- 3 # length of weight vector
constr <- mergeConstraints(
	ordinalConstraint(n, 1, 2),
	ordinalConstraint(n, 1, 3))
transform <- simplex.createTransform(n)
constr <- simplex.createConstraints(transform, constr)
seedPoint <- createSeedPoint(constr, homogeneous=TRUE)
N <- 1000
w <- har(seedPoint, constr, N=N * (n-1)^3, thin=(n-1)^3, homogeneous=TRUE, transform=transform)$samples
stopifnot(all(w[,1] >= w[,2]) && all(w[,1] >= w[,3]))
}
