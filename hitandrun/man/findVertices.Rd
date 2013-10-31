\name{findVertices}
\alias{findVertices}
\title{
Find vertices of the polytope 
}
\description{
Find the vertices of a polytope given by a set of linear constraints.
}
\usage{
findVertices(constr, homogeneous=FALSE)
}
\arguments{
	\item{constr}{Constraint definition}
	\item{homogeneous}{Whether constr is given in homogeneous coordinate representation}
}
\details{
See \code{\link{har}} for a description of the constraint definition and the homogeneous coordinate representation.

Uses the Avis-Fukuda pivoting algorithm to enumerate the vertices of the polytope.
}
\value{
A matrix, in which each row is a vertex of the polytope.
}
\author{
Gert van Valkenhoef
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{har}}

\code{\link{findExtremePoints}}
\code{\link{findInteriorPoint}}

\code{\link[rcdd]{scdd}}
}
\examples{
# constraints: x_1 >= 0, x_2 >= 0, x_1 + x_2 <= 1
A <- rbind(c(-1, 0), c(0, -1), c(1, 1))
b <- c(0, 0, 1)
d <- c("<=", "<=", "<=")
constr <- list(constr=A, rhs=b, dir=d)

findVertices(constr, homogeneous=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{seed point}
