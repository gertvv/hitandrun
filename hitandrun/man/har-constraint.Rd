\name{harConstraints}
\alias{harConstraints}
\alias{simplexConstraints}
\alias{lowerBoundConstraint}
\alias{upperBoundConstraint}
\alias{lowerRatioConstraint}
\alias{upperRatioConstraint}
\alias{exactRatioConstraint}
\alias{ordinalConstraint}
\alias{mergeConstraints}
\title{
Constraint formulation utility functions
}
\description{
These utility functions generate linear constraints
}
\usage{
simplexConstraints(n)
lowerBoundConstraint(n, i, x)
upperBoundConstraint(n, i, x)
lowerRatioConstraint(n, i, j, x)
upperRatioConstraint(n, i, j, x)
exactRatioConstraint(n, i, j, x)
ordinalConstraint(n, i, j)
mergeConstraints(...)
}
\arguments{
	\item{n}{Number of dimensions (vector components)}
	\item{i}{Index of first component}
	\item{j}{Index of second component}
	\item{x}{Scalar bound}
	\item{...}{Constraint definitions, or a single list of constraint definitions}
}
\details{
See \code{\link{har}} for a description of the constraint format.

\code{simplexConstraints} encodes the n-simplex: \eqn{\forall_k w_{k} \geq 0}{all k: w_k >= 0} and \eqn{\sum_k w_k = 1}{sum_k(w_k) = 1}

\code{lowerBoundConstraint} encodes \eqn{w_{i} \geq x}{w_i >= x}

\code{upperBoundConstraint} encodes \eqn{w_{i} \leq x}{w_i <= x}

\code{lowerRatioConstraint} encodes \eqn{w_{i} / w_{j} \geq x}{w_i / w_j >= x}

\code{upperRatioConstraint} encodes \eqn{w_{i} / w_{j} \leq x}{w_i / w_j <= x}

\code{exactRatioConstraint} encodes \eqn{w_{i} / w_{j} = x}{w_i / w_j = x}

\code{ordinalConstraint} encodes \eqn{w_{i} \geq w_{j}}{w_i >= w_j}

\code{mergeConstraints} merges the constraints it is given. Alternatively, the function takes a single list of constraint definitions which are to be merged.
}
\value{
A constraint definition (concatenation of the given constraint definitions).
}
\author{
Gert van Valkenhoef
}
\seealso{
\code{\link{eliminateRedundant}}
\code{\link{hitandrun}}
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
