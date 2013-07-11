\name{createBoundBox}
\alias{createBoundBox}
\title{
Calculate a bounding box
}
\description{
Calculate a bounding box around a polytope given by a set of linear constraints.
}
\usage{
createBoundBox(constr, homogeneous=FALSE)
}
\arguments{
	\item{constr}{Constraint definition}
	\item{homogeneous}{Whether constr is given in homogeneous coordinate representation}
}
\details{
See \code{\link{har}} for a description of the constraint definition and the homogeneous coordinate representation.

This function uses \code{\link{findExtremePoints}} to find extreme points along each dimension.
}
\value{
\item{lb}{Lower bound for each dimension (not including homogeneous coordinate).}
\item{ub}{Upper bound for each dimension (not including homogeneous coordinate).}
}
\author{
Gert van Valkenhoef
}
\seealso{
\code{\link{har}}

\code{\link{findExtremePoints}}
}
\examples{
# constraints: x_1 >= 0, x_2 >= 0, x_1 + x_2 <= 1
A <- rbind(c(-1, 0), c(0, -1), c(1, 1))
b <- c(0, 0, 1)
d <- c("<=", "<=", "<=")
constr <- list(constr=A, rhs=b, dir=d)

bb <- createBoundBox(constr)
stopifnot(bb$lb == c(0.0, 0.0))
stopifnot(bb$ub == c(1.0, 1.0))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{bounding box}
