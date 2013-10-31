\name{simplex.createConstraints}
\alias{simplex.createConstraints}
\title{
Create constraints that define the (n-1)-simplex
}
\description{
This function takes a transformation matrix from the plane coincident with the (n-1) simplex and (optionally) additional constraints defined in n-dimensional space, and generates a set of constraints defining the simplex and (optionally) the additional constraints in the (n-1)-dimensional homogeneous coordinate system.
}
\usage{
simplex.createConstraints(transform, userConstr=NULL)
}
\arguments{
	\item{transform}{Transformation matrix}
	\item{userConstr}{Additional constraints}
}
\details{
The transformation of the constraint matrix to (n-1)-dimensional homogeneous coordinates is a necessary preprocessing step for applying "Hit and Run" to subsets of the simplex defined by \code{userConstr}.
}
\value{
  A set of constraints in the (n-1)-dimensional homogeneous coordinate system.
}
\author{
Gert van Valkenhoef
}
\seealso{
\code{\link{simplex.createTransform}}
\code{\link{har}}
\code{\link{harConstraints}}
}
\examples{
n <- 3
userConstr <- mergeConstraints(
	ordinalConstraint(3, 1, 2), ordinalConstraint(3, 2, 3))

transform <- simplex.createTransform(n)
constr <- simplex.createConstraints(transform, userConstr)
seedPoint <- createSeedPoint(constr, homogeneous=TRUE)

N <- 10000
samples <- har(seedPoint, constr, N, 1, homogeneous=TRUE, transform=transform)$samples

# Check dimension
stopifnot(dim(samples) == c(N, n))

# Check that w_i >= w_i+1
stopifnot(sapply(1:(n-1), function(i) {
	all(samples[,i]>=samples[,i+1])
}))

# Check that w_i >= 0
stopifnot(samples >= 0)

# Check that sum_i w_i = 1
E <- 1E-12
stopifnot(apply(samples, 1, sum) > 1 - E)
stopifnot(apply(samples, 1, sum) < 1 + E)
}
\keyword{simplex}
\keyword{constraint}
\keyword{constraint}
