# create the constraint that w1 >= w2 in an n-dim space
ordinalConstraint <- function(n, w1, w2) {
	a <- rep(0, n + 1)
	a[w1] <- -1
	a[w2] <- 1
	t(a)
}

# create the constraint that w1 <= x
upperBoundConstraint <- function(n, w1, x) {
	a <- rep(0, n + 1)
	a[w1] <- 1
	a[n + 1] <- x
	t(a)
}

# create the constraint that w1 >= x
lowerBoundConstraint <- function(n, w1, x) {
	a <- rep(0, n + 1)
	a[w1] <- -1
	a[n + 1] <- -x
	t(a)
}

# create the constraint that w1/w2 <= x
upperRatioConstraint <- function(n, w1, w2, x) {
	a <- rep(0, n + 1)
	a[w1] <- 1
	a[w2] <- -x
	t(a)
}

# create the constraint that x <= w1/w2
lowerRatioConstraint <- function(n, w1, w2, x) {
	a <- rep(0, n + 1)
	a[w1] <- -1
	a[w2] <- x
	t(a)
}
