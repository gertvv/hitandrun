\name{solution.basis}
\alias{solution.basis}
\title{
Calculate the basis for the solution space of a system of linear equations
}
\description{
Given a set of linear equality constraints, determine a translation and a basis for its solution space.
}
\usage{
solution.basis(constr)
}
\arguments{
	\item{constr}{Linear equality constraints}
}
\details{
For a system of linear equations, \eqn{Ax = b}{Ax = b}, the solution space is given by
\deqn{x = A^\dagger b + (I - A^\dagger A) y}{x = A'b + (I - A' A)}
where \eqn{A^\dagger}{A'} is the Moore-Penrose pseudoinverse of \eqn{A}{A}.
The QR decomposition of \eqn{I - A^\dagger A}{I - A'A} enables us to determine the dimension of the solution space and derive a basis for that space.
}
\value{
A list, consisting of 
\item{translate}{A point in the solution space}
\item{basis}{A basis rooted in that point}
}
\author{
Gert van Valkenhoef
}
\seealso{
\code{\link{createTransform}}
}
\examples{
# A 3-dimensional original space
n <- 3

# x_1 + x_2 + x_3 = 1
eq.constr <- list(constr = t(rep(1, n)), dir = '=', rhs = 1)
basis <- solution.basis(eq.constr)
stopifnot(ncol(basis$basis) == 2) # Dimension reduced to 2
y <- rbind(rnorm(100, 0, 100), rnorm(100, 0, 100))
x <- basis$basis \%*\% y + basis$translate
stopifnot(all.equal(apply(x, 2, sum), rep(1, 100)))

# 2 x_2 = x_1; 2 x_3 = x_2
eq.constr <- mergeConstraints(
  eq.constr,
  list(constr = c(-1, 2, 0), dir = '=', rhs = 0),
  list(constr = c(0, -1, 2), dir = '=', rhs = 0))
basis <- solution.basis(eq.constr)
stopifnot(ncol(basis$basis) == 0) # Dimension reduced to 0
stopifnot(all.equal(basis$translate, c(4/7, 2/7, 1/7)))
}
\keyword{transform}
