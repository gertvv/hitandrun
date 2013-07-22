# create the constraint that w_i1 >= w_i2
ordinalConstraint <- function(n, i1, i2) {
  a <- rep(0, n)
  a[i1] <- -1
  a[i2] <- 1
  list(constr=t(a), rhs=c(0), dir=c("<="))
}

# create the constraint that w_i1 <= x
upperBoundConstraint <- function(n, i1, x) {
  a <- rep(0, n)
  a[i1] <- 1
  a[n + 1] <- x
  list(constr=t(a), rhs=c(0), dir=c("<="))
}

# create the constraint that w_i1 >= x
lowerBoundConstraint <- function(n, i1, x) {
  a <- rep(0, n)
  a[i1] <- -1
  a[n + 1] <- -x
  list(constr=t(a), rhs=c(0), dir=c("<="))
}

# create the constraint that w_i1/w_i2 <= x
upperRatioConstraint <- function(n, i1, i2, x) {
  a <- rep(0, n)
  a[i1] <- 1
  a[i2] <- -x
  list(constr=t(a), rhs=c(0), dir=c("<="))
}

# create the constraint that x <= w_i1/w_i2
lowerRatioConstraint <- function(n, i1, i2, x) {
  a <- rep(0, n)
  a[i1] <- -1
  a[i2] <- x
  list(constr=t(a), rhs=c(0), dir=c("<="))
}

# merge a list of constraints
mergeConstraints <- function(...) {
  lst <- list(...)
  if (length(lst) == 1 && is.list(lst[[1]])) {
    lst <- lst[[1]]
  }
  list(
    constr=do.call("rbind", lapply(lst, function(c) { c$constr })),
    rhs=unlist(lapply(lst, function(c) { c$rhs })),
    dir=unlist(lapply(lst, function(c) { c$dir }))
  )
}
