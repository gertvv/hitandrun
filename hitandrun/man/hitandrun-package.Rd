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
}
