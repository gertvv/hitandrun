\name{eliminateRedundant}
\alias{eliminateRedundant}
\title{
Eliminate redundant linear constraints
}
\description{
Given a set of linear constraints, gives a subset of these constraints that are non-redundant.
}
\usage{
eliminateRedundant(constr)
}
\arguments{
	\item{constr}{Constraints}
}
\details{
  If no constraints are redundant, returns the same set of constraints.
}
\value{
  A set of non-redundant constraints.
}
\author{
  Gert van Valkenhoef, Tommi Tervonen
}
\seealso{
\code{\link{harConstraints}}
}
\examples{
constr <- list(
  constr = rbind(
    c(-1  ,  0),
    c( 0  , -1),
    c( 1  ,  1),
    c( 0.5, -1)),
  dir = c('<=', '<=', '=', '<='),
  rhs = c(0, 0, 1, 0))

constr <- eliminateRedundant(constr)

stopifnot(nrow(constr$constr) == 3) # eliminates one constraint
}
\keyword{constraint}
