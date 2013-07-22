library('hitandrun')

## Integration test for bug #2: https://github.com/gertvv/hitandrun/issues/2

n <- 4
ord <- c(4,2,1,3)
pairs <- cbind(ord[1:(n-1)], ord[2:n])
constr <- mergeConstraints(
  apply(pairs, 1, function(pair) {
    ordinalConstraint(n, pair[1], pair[2])
  }))
transform <- simplex.createTransform(n)
constr <- simplex.createConstraints(transform, constr)
seedPoint <- createSeedPoint(constr, homogeneous=TRUE)

print(seedPoint %*% t(transform))
stopifnot(all((constr$constr %*% seedPoint) <= constr$rhs))

## Integration test for bug #1: https://github.com/gertvv/hitandrun/issues/1
# "Generated seed point may be outside polytope"

# 1 > 4 > 5 > 6 > {2, 3}
constr <- structure(list(constr = structure(c(-1, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 0, 1, 1, -1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1,
  -1, -1), .Dim = 5:6), rhs = c(0, 0, 0, 0, 0), dir = c("<=", "<=",
  "<=", "<=", "<=")), .Names = c("constr", "rhs", "dir"))
n <- ncol(constr$constr)

transform <- simplex.createTransform(n)
constr <- simplex.createConstraints(transform, constr)
seedPoint <- createSeedPoint(constr, homogeneous=TRUE)

print(seedPoint %*% t(transform))
stopifnot(all((constr$constr %*% seedPoint) <= constr$rhs))

