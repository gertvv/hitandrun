## Integration test for bug #1: https://github.com/gertvv/hitandrun/issues/1
# "Generated seed point may be outside polytope"

library('hitandrun')

constr <- structure(list(constr = structure(c(-1, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 0, 1, 1, -1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1,
  -1, -1), .Dim = 5:6), rhs = c(0, 0, 0, 0, 0), dir = c("<=", "<=",
  "<=", "<=", "<=")), .Names = c("constr", "rhs", "dir"))
n <- ncol(constr$constr)

transform <- simplex.createTransform(n)
constr <- simplex.createConstraints(transform, constr)
seedPoint <- createSeedPoint(constr, homogeneous=TRUE)

stopifnot(all((constr$constr %*% seedPoint) <= constr$rhs))
