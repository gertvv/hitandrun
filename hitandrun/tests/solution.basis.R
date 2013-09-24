library(hitandrun)

# A 3-dimensional original space
n <- 3

# x_1 + x_2 + x_3 = 1
eq.constr <- list(
  constr = matrix(rep(1, n), nrow=1, ncol=n),
  dir = '=',
  rhs = 1)

basis <- solution.basis(eq.constr)

stopifnot(all.equal(sum(basis$translate), 1))
stopifnot(nrow(basis$basis) == 3)
stopifnot(ncol(basis$basis) == 2) # Dimension reduced to 2

y <- rbind(rnorm(100, 0, 100), rnorm(100, 0, 100))
x <- basis$basis %*% y + basis$translate
stopifnot(all.equal(apply(x, 2, sum), rep(1, 100)))

# 2 x_2 = x_1
eq.constr <- mergeConstraints(
  eq.constr,
  list(constr = c(-1, 2, 0), dir='=', rhs=0))

basis <- solution.basis(eq.constr)
stopifnot(all.equal(sum(basis$translate), 1))
stopifnot(all.equal(basis$translate[1], 2 * basis$translate[2]))
stopifnot(nrow(basis$basis) == 3)
stopifnot(ncol(basis$basis) == 1) # Dimension reduced to 1

y <- y[1,,drop=FALSE]
x <- basis$basis %*% y + basis$translate
stopifnot(all.equal(apply(x, 2, sum), rep(1, 100)))
stopifnot(all.equal(x[1,], 2 * x[2,]))

# 2 x_3 = x_2
eq.constr <- mergeConstraints(
  eq.constr,
  list(constr = c(0, -1, 2), dir='=', rhs=0))

basis <- solution.basis(eq.constr)
stopifnot(all.equal(sum(basis$translate), 1))
stopifnot(all.equal(basis$translate[1], 2 * basis$translate[2]))
stopifnot(all.equal(basis$translate[2], 2 * basis$translate[3]))
stopifnot(nrow(basis$basis) == 3)
stopifnot(ncol(basis$basis) == 0) # Dimension reduced to 0
