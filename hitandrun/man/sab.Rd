\name{sab}
\alias{sab}
\title{
"Shake and Bake" sampler
}
\description{
The "Shake and Bake" method generates a Markov Chain whose stable state converges on the uniform distribution over the boundary of a convex polytope defined by a set of linear constraints.
}
\usage{
sab(x0, i0, constr, N, thin=1, homogeneous=FALSE, transform=NULL)
}
\arguments{
	\item{x0}{Starting point (must be in the polytope)}
	\item{i0}{Index of the closest face to the starting point}
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
See \code{\link{shakeandbake}} for a "Shake and Bake" sampler that also supports equality constraints.
The constraints define the polytope as usual for linear programming: \eqn{Ax \leq b}{Ax <= b}.
In particular, it must be true that \eqn{A x_0 \leq b}{Ax0 <= b}.
Points are generated from the boundary of the polytope (where equality holds for one of the constraints), using the "running" shake and bake sampler, which samples the direction vector so that every move point is accepted (Boender et al. 1991).
}
\value{
A list, containing:
	\item{samples}{A matrix containing the generated samples as rows.}
	\item{faces}{A vector containing the indices of the faces on which the samples lie.}
	\item{xN}{The last generated sample, untransformed. Can be used as the starting point for a continuation of the chain.}
	\item{iN}{Face on which the last generated sample lies.}
}
\references{
	Boender, C. G. E., Caron, R. J., McDonald, J. F., Rinnooy Kan, A. H. G., Romeijn, H. E., Smith, R. L., Telgen, J., and Vorst, A. C. F. (1991) "Shake-and-Bake Algorithms for Generating Uniform Points on the Boundary of Bounded Polyhedra". \emph{Operations Research} 39(6):945-954. \href{http://dx.doi.org/10.1287/opre.39.6.945}{doi:10.1287/opre.39.6.945}
}
\author{
Gert van Valkenhoef
}
\note{
"Shake and Bake" is a Markov Chain Monte Carlo (MCMC) method, so generated samples form a correlated time series.
}
\seealso{
\code{\link{harConstraints}}
\code{\link{shakeandbake}}
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
result <- sab(x0, 1, constr, 1E4)
samples <- result$samples

# Check dimension of result
stopifnot(dim(samples) == c(1E4, 2))

# Check that x_i >= 0
stopifnot(samples >= -1E-15)

# Check that x_1 + x_2 <= 1
stopifnot(samples[,1] + samples[,2] <= 1 + 1E-15)

# check that the results lie on the faces
faces <- result$faces
stopifnot(all.equal(samples[faces==1,1], rep(0, sum(faces==1))))
stopifnot(all.equal(samples[faces==2,2], rep(0, sum(faces==2))))
stopifnot(all.equal(samples[faces==3,1] + samples[faces==3,2], rep(1, sum(faces==3))))

\dontrun{plot(samples)}

}
\keyword{shake-and-bake}
\keyword{uniform sampling}
