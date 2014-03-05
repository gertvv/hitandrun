\name{eliminateRedundant}
\alias{eliminateRedundant}
\title{
  Eliminate redundant linear constraints
}
\description{
Given a set of linear constraints, gives a subset of these constraints
that are non-redundant, i.e. that are not required for defining the
convex polytope.
}
\usage{
eliminateRedundant(constr)
}
\arguments{
	\item{constr}{A structure of HAR linear constraints}
}
\details{
  If no constraints are redundant, returns the same set of constraints.
}
\value{
  A set of non-redundant constraints.
}
\author{
  Gert van Valkenhoef and Tommi Tervonen <tommi@smaa.fi>
}
\seealso{
  \code{\link{har}}
}
\examples{

constr <- list(
  constr = matrix(c(-1, 0, 0, -1, 1, 1, 0.5, -1), ncol=2, byrow=TRUE),
  dir = c('<=', '<=', '=', '<='),
  rhs = c(0, 0, 1, 0))

constr <- eliminateRedundant(constr)

stopifnot(nrow(constr$constr) == 3) # eliminates one constraint
}

\keyword{constraint}
