\name{hitandrun-package}
\alias{hitandrun-package}
\docType{package}
\title{
"Hit and Run" sampling
}
\description{
This package provides a "Hit and Run" sampler that generates a Markov chain whose stable state converges on the uniform distribution over a convex polytope.
The polytope is given by a set of inequality constraints in standard linear programming form (\eqn{Ax \leq b}{Ax <= b}) and optionally a set of equality constraints.
In addition, there is a "Shake and Bake" sampler to generate points from the boundary of such a shape.

Utilities are provided for sampling from subsets of the unit simplex (i.e. random variates that can be interpreted as weights satisfying certain constraints) and for specifying common constraints.
}
\details{
\code{\link{hitandrun}} and \code{\link{shakeandbake}} now provide the most general interface for sampling from spaces defined by arbitrary linear equality and inequality constraints. The functions described in the following provide lower level functionality on which it is built.

\code{\link{har}} is the core "Hit and Run" sampler, \code{\link{sab}} is the core "Shake and Bake" sampler, \code{\link{bbReject}} is the bounding box rejection sampler, and \code{\link{simplex.sample}} samples uniformly from the unit simplex.

See \code{\link{simplex.createTransform}} and \code{\link{simplex.createConstraints}} for sampling from subsets of the unit simplex.
Utilities to specify common constraints are described in \code{\link{harConstraints}}.

When the sampling space is restricted by different linear equality constraints, use \code{\link{solution.basis}}, \code{\link{createTransform}}, and \code{\link{transformConstraints}}.
This is a generalization of the methods for sampling from the simplex.
}
\author{
Maintainer: Gert van Valkenhoef <g.h.m.van.valkenhoef@rug.nl>
}
\references{
    Tervonen, T., van Valkenhoef, G., Basturk, N., and Postmus, D. (2012) "Hit-And-Run enables efficient weight generation for simulation-based multiple criteria decision analysis". \emph{European Journal of Operational Research} 224(3) 552-559. \href{http://dx.doi.org/10.1016/j.ejor.2012.08.026}{doi:10.1016/j.ejor.2012.08.026}
    van Valkenhoef, G., Tervonen, T., and Postmus, D. (2014) "Notes on 'Hit-And-Run enables efficient weight generation for simulation-based multiple criteria decision analysis'". \emph{European Journal of Operational Research} (in press). \href{http://dx.doi.org/10.1016/j.ejor.2014.06.036}{doi:10.1016/j.ejor.2014.06.036}
}
\note{
"Hit and Run" is a Markov Chain Monte Carlo (MCMC) method, so generated samples form a correlated time series. To get a uniform sample, you need \eqn{O^*(n^3)}{O*(n^3)} samples, where n is the dimension of the sampling space.
}
\keyword{hit-and-run}
\seealso{
\code{\link{hitandrun}}
\code{\link{har}}

\code{\link{bbReject}}
\code{\link{simplex.sample}}
\code{\link{hypersphere.sample}}

\code{\link{solution.basis}}
\code{\link{createTransform}}
\code{\link{transformConstraints}}

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
w <- har(seedPoint, constr, N=N * (n-1)^3, thin=(n-1)^3,
  homogeneous=TRUE, transform=transform)$samples
stopifnot(all(w[,1] >= w[,2]) && all(w[,1] >= w[,3]))
}
