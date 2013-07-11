\name{createSeedPoint}
\alias{createSeedPoint}
\title{
Generate a seed point
}
\description{
Generate a seed point inside a polytope given by a set of linear constraints.
}
\usage{
createSeedPoint(constr, homogeneous=FALSE, randomize=FALSE, method="extremes")
}
\arguments{
	\item{constr}{Constraint definition}
	\item{homogeneous}{Whether constr is given in homogeneous coordinate representation}
	\item{randomize}{If TRUE, randomly weight the generating points}
	\item{method}{How to obtain generating points: "extremes" for extreme points along each dimension, "vertices" for the vertices of the polytope}
}
\details{
See \code{\link{har}} for a description of the constraint definition and the homogeneous coordinate representation.

This method generates a point within the polytope defined by the given constraints by calculating the weighted arithmetic mean of a set of generating points.
If `randomize' is set, the weights are randomly generated, otherwise they are all equal.

The "extremes" method uses \code{\link{findExtremePoints}} to find extreme points along each dimension, and the "vertices" method uses \code{\link{findVertices}} to find the vertices of the polytope.
The "vertices" method is computationally more expensive.
}
\value{
A coordinate vector in the appropriate coordinate system.
}
\author{
Gert van Valkenhoef
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{har}}

\code{\link{findExtremePoints}}
\code{\link{findVertices}}
}
\examples{
# constraints: x_1 >= 0, x_2 >= 0, x_1 + x_2 <= 1
A <- rbind(c(-1, 0), c(0, -1), c(1, 1))
b <- c(0, 0, 1)
d <- c("<=", "<=", "<=")
constr <- list(constr=A, rhs=b, dir=d)

x0 <- createSeedPoint(constr)
stopifnot(x0 >= 0)
stopifnot(sum(x0) <= 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{seed point}
