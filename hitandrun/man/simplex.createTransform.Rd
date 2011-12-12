\name{simplex.createTransform}
\alias{simplex.createTransform}
\title{
Transform points on an (n-1)-simplex to n-dimensional space
}
\description{
This function generates a matrix to transform points in an (n-1) dimensional homogeneous coordinate representation of the (n-1) simplex to n-dimensional Cartesian coordinates.
}
\usage{
simplex.createTransform(n)
}
\arguments{
	\item{n}{Dimension of the space}
}
\details{
}
\value{
  A transformation matrix.
}
\references{
}
\author{
Gert van Valkenhoef
}
\note{
%%  ~~further notes~~
}
\seealso{
\code{\link{simplex.createConstraints}}
\code{\link{har}}
}
\examples{
transform <- simplex.createTransform(3)
# The origin in (n-1)-dimensional space should lie in the simplex when
# transformed to n-dimensional space
x <- transform \%*\% c(0, 0, 1)

# Check that x_i >= 0
stopifnot(x >= 0)
# Check that \sum_i x_i = 1
E <- 1E-12
stopifnot(sum(x) > 1 - E)
stopifnot(sum(x) < 1 + E)
}
\keyword{simplex}
