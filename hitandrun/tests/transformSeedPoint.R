# In versions <= 0.4-1, no transformation was applied to the supplied starting point by
# har.init() (and therefore hitandrun()).

library(hitandrun)

constr <- list(
  constr=rbind(rep(1, 3), -diag(3)),
  rhs=c(1, rep(0, 3)),
  dir=c("=", rep("<=", 3)))
x0 <- c(1/2,1/4,1/4)

state <- har.init(constr, x0=x0)

stopifnot(all.equal(as.vector(state$transform %*% state$x0), x0))

constr$dir <- rep("<=", 4)
state <- har.init(constr, x0=x0)
print(state$x0)
stopifnot(all.equal(as.vector(state$x0), c(x0, 1)))
stopifnot(all.equal(as.vector(state$transform %*% state$x0), x0))
