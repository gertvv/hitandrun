\name{createTransform}
\alias{createTransform}
\title{
Calculate the basis for the solution space of a system of linear equations
}
\description{
This function takes a basis, consisting of an \eqn{n \times m}{n*m} change of basis matrix and an \eqn{n}{n}-vector representing the origin of the \eqn{m}{m}-space, and generates a matrix to transform points in the \eqn{m}{m}-space, given in homogeneous coordinates, to the \eqn{n}-space.

The inverse transform can also be generated, and conversion can be to homogeneous coordinates instead of Cartesion ones.
}
\usage{
createTransform(basis, inverse=FALSE, keepHomogeneous=inverse)
}
\arguments{
	\item{basis}{Basis (and origin) for the \eqn{m}{m}-space (see \code{\link{solution.basis}})}
	\item{inverse}{TRUE to convert from \eqn{n}{n}-space coordinates to \eqn{m}{m}-space coordinates}
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
\code{\link{solution.basis}}
}
\keyword{transform}
