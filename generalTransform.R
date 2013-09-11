library(hitandrun)

runHar <- function(eq.constr, ineq.constr) {
  n <- ncol(eq.constr$constr)
  stopifnot(ncol(ineq.constr$constr) == n)

  basis <- solution.basis(eq.constr)
  transform <- createTransform(basis)
  constr <- simplex.createConstraints(transform, ineq.constr)
  x0 <- createSeedPoint(constr, homogeneous=TRUE)
  print(x0)
  if (ncol(basis$basis) > 0) {
    har(x0, constr, 1E4, transform=transform, homogeneous=TRUE)$samples
  } else {
    matrix(rep(basis$translate, each=1E4), nrow=1E4)
  }
}

n <- 3

eq.constr <- list(
  constr = matrix(rep(1, n), nrow=1, ncol=n),
  dir = '=',
  rhs = 1)

ineq.constr <- mergeConstraints(
  ordinalConstraint(3, 1, 2),
  ordinalConstraint(3, 2, 3))

w <- runHar(eq.constr, ineq.constr)
stopifnot(all(w[,1] >= w[,2]) && all(w[,1] >= w[,3]))
stopifnot(all(w >= 0))
stopifnot(abs(w[,1] + w[,2] + w[,3] - 1) <= 1E-15)

eq.constr <- mergeConstraints(
  eq.constr,
  list(constr = c(-1, 2, 0), dir='=', rhs=0))

w <- runHar(eq.constr, ineq.constr)
stopifnot(all(w[,1] >= w[,2]) && all(w[,1] >= w[,3]))
stopifnot(all(w >= 0))
stopifnot(abs(w[,1] + w[,2] + w[,3] - 1) <= 1E-15)
stopifnot(abs(w[,1] - 2 * w[,2]) <= 1E-10)

eq.constr <- mergeConstraints(
  eq.constr,
  list(constr = c(0, -1, 2), dir='=', rhs=0))

w <- runHar(eq.constr, ineq.constr)
