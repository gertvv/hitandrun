\name{createSeedPoint}
\alias{createSeedPoint}
\title{
Generate a seed point
}
\description{
Generate a seed point inside a polytope given by a set of linear constraints.
}
\usage{
createSeedPoint(constr, homogeneous=FALSE, randomize=FALSE, method="slacklp")
}
\arguments{
	\item{constr}{Constraint definition}
	\item{homogeneous}{Whether constr is given in homogeneous coordinate representation}
	\item{randomize}{If TRUE, randomize the starting point}
  \item{method}{How to obtain the starting point: "slacklp" for a linear program that maximizes the minimum slack, or "vertices" for a weighted average of the vertices of the polytope}
}
\details{
See \code{\link{har}} for a description of the constraint definition and the homogeneous coordinate representation.


\itemize{
\item{The "slacklp" method solves a linear program that maximizes the minimum slack on the inequality constraints.
When randomized, the slack on each constraint is randomly rescaled before maximization.}
\item{The "vertices" method enumerates all vertices of the polytope and then calculates the weighted arithmetic mean of this set of points. If `randomize' is set, the weights are randomly generated, otherwise they are all equal and the generated point is the centroid of the polytope.}
}
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
