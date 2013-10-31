\name{findExtremePoints}
\alias{findExtremePoints}
\title{
Find extreme points
}
\description{
Find extreme points of a polytope given by a set of linear constraints along each dimension.
}
\usage{
findExtremePoints(constr, homogeneous=FALSE)
}
\arguments{
	\item{constr}{Constraint definition}
	\item{homogeneous}{Whether constr is given in homogeneous coordinate representation}
}
\details{
See \code{\link{har}} for a description of the constraint definition and the homogeneous coordinate representation.

For n-dimensional coordinate vectors, solves 2n LPs to find the extreme points along each dimension.
}
\value{
A matrix, in which each row is a coordinate vector in the appropriate coordinate system.
}
\author{
Gert van Valkenhoef
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{har}}

\code{\link{findInteriorPoint}}
\code{\link{findVertices}}

\code{\link[rcdd]{lpcdd}}
}
\examples{
# constraints: x_1 >= 0, x_2 >= 0, x_1 + x_2 <= 1
A <- rbind(c(-1, 0), c(0, -1), c(1, 1))
b <- c(0, 0, 1)
d <- c("<=", "<=", "<=")
constr <- list(constr=A, rhs=b, dir=d)

findExtremePoints(constr, homogeneous=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{seed point}
