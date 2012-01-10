\name{harConstraints}
\alias{harConstraints}
\alias{lowerBoundConstraint}
\alias{upperBoundConstraint}
\alias{lowerRatioConstraint}
\alias{upperRatioConstraint}
\alias{ordinalConstraint}
\alias{mergeConstraints}
\title{
Constraint formulation utility functions
}
\description{
These utility functions generate constraints
}
\usage{
lowerBoundConstraint(n, i1, x)
upperBoundConstraint(n, i1, x)
lowerRatioConstraint(n, i1, i2, x)
upperRatioConstraint(n, i1, i2, x)
ordinalConstraint(n, i1, i2)
mergeConstraints(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{n}{Number of dimensions (vector components)}
	\item{i1}{Index of first component}
	\item{i2}{Index of second component}
	\item{x}{Scalar bound}
	\item{...}{Constraint definitions, or a single list of constraint definitions}
}
\details{
See \code{\link{har}} for a description of the constraint format.

\code{lowerBoundConstraint} encodes w_i1 >= x

\code{upperBoundConstraint} encodes w_i1 <= x

\code{lowerRatioConstraint} encodes w_i1 / w_i2 >= x

\code{upperRatioConstraint} encodes w_i1 / w_i2 <= x

\code{ordinalConstraint} encodes w_i1 >= w_i2

\code{mergeConstraints} merges the constraints it is given. Alternatively, the function takes a single list of constraint definitions which are to be merged.
}
\value{
A constraint definition (concatenation of the given constraint definitions).
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Gert van Valkenhoef
}
\note{
%%  ~~further notes~~
}
\seealso{
\code{\link{har}}
}
\examples{
# create an ordinal constraint
c1 <- ordinalConstraint(2, 1, 2)
stopifnot(c1$constr == c(-1, 1))
stopifnot(c1$rhs == c(0))
stopifnot(c1$dir == c("<="))

# create our own constraints
c2 <- list(constr=t(c(-1, 0)), rhs=c(0), dir=c("<="))
c3 <- list(constr=t(c(1, 1)), rhs=c(1), dir=c("<="))

# merge the constraints into a single definition
c <- mergeConstraints(c1, c2, c3)
stopifnot(c$constr == rbind(c(-1, 1), c(-1, 0), c(1, 1)))
stopifnot(c$rhs == c(0, 0, 1))
stopifnot(c$dir == c("<=", "<=", "<="))

# test the alternative (list) method
l <- mergeConstraints(list(c1, c2, c3))
stopifnot(c$constr == l$constr)
stopifnot(c$rhs == l$rhs)
stopifnot(c$dir == l$dir)

# test iteratively merging
l <- mergeConstraints(mergeConstraints(c1, c2), c3)
stopifnot(c$constr == l$constr)
stopifnot(c$rhs == l$rhs)
stopifnot(c$dir == l$dir)
}
\keyword{hit-and-run}
\keyword{constraint}
