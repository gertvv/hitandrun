\name{hitandrun}
\alias{hitandrun}
\title{
"Hit and Run" sampler
}
\description{
The "Hit and Run" method generates a Markov Chain whose stable state converges on the uniform distribution over a convex polytope defined by a set of linear inequality constraints. \code{hitandrun()} further uses the Moore-Penrose pseudo-inverse to eliminate an arbitrary set of linear equality constraints before applying the "Hit and Run" sampler.
}
\usage{
hitandrun(constr, n.samples=1E4,
    thin.fn = function(n) { ceiling(log(n + 1)/4 * n^3) }, thin = NULL,
    x0 = NULL)
}
\arguments{
  \item{constr}{Linear constraints that define the sampling space (see details)}
	\item{n.samples}{The desired number of samples to return. The sampler is run for \code{n.samples * thin} iterations}
	\item{thin.fn}{Function that specifies a thinning factor depending on the dimension of the sampling space after equality constraints have been eliminated. Will only be invoked if \code{thin} is \code{NULL}}
  \item{thin}{The thinning factor}
  \item{x0}{Seed point for the Markov Chain. Can be used to continue sampling if convergence is inadequate (see value)}
}
\details{
The constraints are given as a list with the elements \code{constr}, \code{dir} and \code{rhs}. \code{dir} is a vector with values \code{'='} or \code{'<='}. \code{constr} is a matrix and \code{rhs} a vector, which encode the standard linear programming constraint froms \eqn{Ax = b} and \eqn{Ax \leq b} (depending on \code{dir}). The lengths of \code{rhs} and \code{dir} must match the number of rows of \code{constr}.

\code{hitandrun()} applies \code{\link{solution.basis}} to generate a basis of the (translated) solution space of the linear constraints (if any). An affine transformation is generated using \code{\link{createTransform}} and applied to the constraints. Then, a seed point satisfying the inequality constraints is generated using \code{\link{createSeedPoint}}. Finally, \code{\link{har}} is used to generate the samples.
}
\value{
A list, containing:
	\item{samples}{A matrix containing the generated samples as rows.}
	\item{xN}{The last generated sample, untransformed. Can be used as the starting point for a continuation of the chain.}
}
\references{
	Smith, R. L. (1984) "Efficient Monte Carlo Procedures for Generating Points Uniformly Distributed over Bounded Regions". \emph{Operations Research} 32(6): 1296-1308. \href{http://dx.doi.org/10.1287/opre.32.6.1296}{doi:10.1287/opre.32.6.1296}
}
\author{
Gert van Valkenhoef
}
\note{
"Hit and Run" is a Markov Chain Monte Carlo (MCMC) method, so generated samples form a correlated time series. To get a uniform sample, you need O*(n^3) samples, where n is the dimension of the sampling space.
}
\seealso{
\code{\link{harConstraints}}
\code{\link{har}}
}
\examples{
# Sample from the 3-simplex with the additional constraint that w_1/w_2 = 2
# Three inequality constraints, two equality constraints
constr <- mergeConstraints(simplexConstraints(3), exactRatioConstraint(3, 1, 2, 2))
result <- hitandrun(constr, n.samples=1000)
samples <- result$samples
stopifnot(dim(samples) == c(1000, 3))
stopifnot(all.equal(apply(samples, 1, sum), rep(1, 1000)))
stopifnot(all.equal(samples[,1]/samples[,2], rep(2, 1000)))

# Sample from the unit rectangle (no equality constraints)
constr <- list(
  constr = rbind(c(1,0), c(0,1), c(-1,0), c(0,-1)),
  dir=rep('<=', 4),
  rhs=c(1, 1, 0, 0))
result <- hitandrun(constr, n.samples=1000)
samples <- result$samples
stopifnot(all(samples >= 0 & samples <= 1))
# Continue sampling from the same chain:
result <- hitandrun(constr, n.samples=1000, x0=result$xN)
}
\keyword{hit-and-run}
\keyword{uniform}
\keyword{sampling}
