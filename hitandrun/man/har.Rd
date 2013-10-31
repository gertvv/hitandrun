\name{har}
\alias{har}
\title{
"Hit and Run" sampler
}
\description{
The "Hit and Run" method generates a Markov Chain whose stable state converges on the uniform distribution over a convex polytope defined by a set of linear constraints.
}
\usage{
har(x0, constr, N, thin=1, homogeneous=FALSE, transform=NULL)
}
\arguments{
	\item{x0}{Starting point (must be in the polytope)}
	\item{constr}{Constraint definition (see details)}
	\item{N}{Number of iterations to run}
	\item{thin}{Thinning factor (keep every 'thin'-th sample)}
	\item{homogeneous}{Whether x0, constr and transform are given in homogeneous coordinate representation (see details)}
	\item{transform}{Transformation matrix to apply to the generated samples (optional)}
}
\details{
The constraints, starting point and transformation matrix can be given in homogeneous coordinate representation (an extra component is added to each vector, equal to 1.0).
This enables affine transformations (such as translation) to be applied to the coordinate vectors by the constraint and transformation matrices.
Be aware that while non-affine (perspective) transformations are also possible, they will not in general preserve uniformity of the generated samples.

Constraints are given as a list(constr=A, rhs=b, dir=d), where d should contain only "<=".
See \code{\link{hitandrun}} for a "Hit and Run" sampler that also supports equality constraints.
The constraints define the polytope as usual for linear programming: \eqn{Ax \leq b}{Ax <= b}.
In particular, it must be true that \eqn{A x_0 \leq b}{Ax0 <= b}.
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
"Hit and Run" is a Markov Chain Monte Carlo (MCMC) method, so generated samples form a correlated time series. To get a uniform sample, you need \eqn{O^*(n^3)}{O*(n^3)} samples, where n is the dimension of the sampling space.
}
\seealso{
\code{\link{harConstraints}}
\code{\link{hitandrun}}
}
\examples{
# constraints: x_1 >= 0, x_2 >= 0, x_1 + x_2 <= 1
A <- rbind(c(-1, 0), c(0, -1), c(1, 1))
b <- c(0, 0, 1)
d <- c("<=", "<=", "<=")
constr <- list(constr=A, rhs=b, dir=d)

# take a point x0 within the polytope
x0 <- c(0.25, 0.25)

# sample 10,000 points
samples <- har(x0, constr, 1E4)$samples

# Check dimension of result
stopifnot(dim(samples) == c(1E4, 2))

# Check that x_i >= 0
stopifnot(samples >= 0)

# Check that x_1 + x_2 <= 1
stopifnot(samples[,1] + samples[,2] <= 1)

\dontrun{plot(samples)}

}
\keyword{hit-and-run}
\keyword{uniform sampling}
