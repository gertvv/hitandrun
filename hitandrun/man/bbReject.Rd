\name{bbReject}
\alias{bbReject}
\title{
Bounding box rejection sampler
}
\description{
Generates uniform random variates over a convex polytope defined by a set of linear constraints by generating uniform variates over a bounding box and rejecting those outside the polytope.
}
\usage{
bbReject(lb, ub, constr, N, homogeneous=FALSE, transform=NULL)
}
\arguments{
	\item{lb}{Lower bound for each dimension (not including homogeneous coordinate)}
	\item{ub}{Upper bound for each dimension (not including homogeneous coordinate)}
	\item{constr}{Constraint definition (see details)}
	\item{N}{Number of samples to generate}
	\item{homogeneous}{Whether constr and transform are given in homogeneous coordinate representation (see details)}
	\item{transform}{Transformation matrix to apply to the generated samples (optional)}
}
\details{
See \code{\link{har}} for a description of the constraint definition and the homogeneous coordinate representation.
}
\value{
A list, containing:
	\item{samples}{A matrix containing the generated samples as rows.}
	\item{rejectionRate}{The mean number of samples rejected for each accepted sample.}
}
\author{
Gert van Valkenhoef
}
\seealso{
\code{\link{createBoundBox}}

\code{\link{harConstraints}}
\code{\link{simplex.createTransform}}
\code{\link{simplex.createConstraints}}
}
\examples{
# constraints: x_1 >= 0, x_2 >= 0, x_1 + x_2 <= 1
A <- rbind(c(-1, 0), c(0, -1), c(1, 1))
b <- c(0, 0, 1)
d <- c("<=", "<=", "<=")
constr <- list(constr=A, rhs=b, dir=d)

# create a bounding box that contains the polytope
lb <- c(0, 0)
ub <- c(1, 1)

# sample 10,000 points
samples <- bbReject(lb, ub, constr, 1E4)$samples

# Check dimension of result
stopifnot(dim(samples) == c(1E4, 2))

# Check that x_i >= 0
stopifnot(samples >= 0)

# Check that x_1 + x_2 <= 1
stopifnot(samples[,1] + samples[,2] <= 1)

\dontrun{plot(samples)}

}
\keyword{bounding box}
\keyword{uniform sampling}
