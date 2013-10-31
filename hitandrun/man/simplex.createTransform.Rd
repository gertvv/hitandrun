\name{simplex.createTransform}
\alias{simplex.createTransform}
\title{
Transform points on an (n-1)-simplex to n-dimensional space
}
\description{
This function generates a matrix to transform points in an (n-1) dimensional homogeneous coordinate representation of the (n-1) simplex to n-dimensional Cartesian coordinates.

The inverse transform can also be generated, and conversion can be to homogeneous coordinates instead of Cartesion ones.
}
\usage{
simplex.createTransform(n, inverse=FALSE, keepHomogeneous=inverse)
}
\arguments{
	\item{n}{Dimension of the space}
	\item{inverse}{TRUE to convert from n-space coordinates to (n-1)-simplex coordinates}
	\item{keepHomogeneous}{TRUE to convert to homogeneous coordinates rather than Cartesian}
}
\details{
Multiply a coordinate vector in homogeneous coordinates by pre-multiplying by the generated matrix (see examples).
}
\value{
  A transformation matrix.
}
\author{
Gert van Valkenhoef
}
\seealso{
\code{\link{simplex.createConstraints}}
\code{\link{har}}
}
\examples{
E <- 1E-12 # Allowed numerical error

# The origin in (n-1)-dimensional space should be the centroid of the simplex
# when transformed to n-dimensional space
transform <- simplex.createTransform(3)
x <- transform \%*\% c(0, 0, 1)
x
stopifnot(abs(x - c(1/3, 1/3, 1/3)) < E)

# The same should hold for the inverse transformation
invTransform <- simplex.createTransform(3, inverse=TRUE)
y <- invTransform \%*\% c(1/3, 1/3, 1/3, 1)
y
stopifnot(abs(y - c(0, 0, 1)) < E)

# Of course, an arbitrary weight vector should transform back to itself
transform <- simplex.createTransform(3, keepHomogeneous=TRUE)
x <- c(0.2, 0.5, 0.3, 1.0)
y <- transform \%*\% invTransform \%*\% x
y
stopifnot(abs(y - x) < E)

# And we can apply the tranform to a matrix:
a <- cbind(x, x, x)
b <- transform \%*\% invTransform \%*\% a
b
stopifnot(abs(b - a) < E)
}
\keyword{simplex}
\keyword{transform}
