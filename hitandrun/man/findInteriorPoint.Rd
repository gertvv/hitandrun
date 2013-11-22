\name{findInteriorPoint}
\alias{findInteriorPoint}
\title{
Find an interior point
}
\description{
Find an interior point of a polytope given by a set of linear constraints along each dimension.
}
\usage{
findInteriorPoint(constr, homogeneous=FALSE, randomize=FALSE)
}
\arguments{
	\item{constr}{Constraint definition}
	\item{homogeneous}{Whether constr is given in homogeneous coordinate representation}
	\item{randomize}{Whether the point should be randomized}
}
\details{
See \code{\link{har}} for a description of the constraint definition and the homogeneous coordinate representation.

Solves a slack-maximizing LP to find an interior point of the polytope defined by the given constraints.
The randomized version randomly scales the slack on each (non-redundant) constraint.
}
\value{
A vector.
}
\author{
Gert van Valkenhoef
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{har}}

\code{\link{findExtremePoints}}
\code{\link{findVertices}}

\code{\link[rcdd]{lpcdd}}
}
\examples{
# constraints: x_1 >= 0, x_2 >= 0, x_1 + x_2 <= 1
A <- rbind(c(-1, 0), c(0, -1), c(1, 1))
b <- c(0, 0, 1)
d <- c("<=", "<=", "<=")
constr <- list(constr=A, rhs=b, dir=d)

findInteriorPoint(constr, homogeneous=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{seed point}
