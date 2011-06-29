findExtremePoints <- function(constr) {
	n <- dim(constr$constr)[2]

	# because lp_solve assumes vars to be non-negative, split each dimension:
	# x_i = y_{2i-1} - y_{2i}
	# thus each dimension is coded as the difference of two variables, give them
	obj <- function(i) {
		obj <- rep(0, 2 * n - 1)
		obj[2 * i - 1] <- 1
		obj[2 * i] <- -1
		obj
	}
	mat <- t(sapply(1:(n-1), obj))
	a <- constr$constr %*% rbind(mat, c(rep(0, 2 * n - 2), 1))

	# for each of the (n-1) dimensions, solve 2 LPs to find the min/max
	findExtreme <- function(dir) {
		function(i) {
			lp(dir, obj(i), a, constr$dir, constr$rhs)$solution
		}
	}
	rawExtremes <- cbind(sapply(1:(n-1), findExtreme("min")), sapply(1:(n-1), findExtreme("max")))
	mat %*% rawExtremes
}

# create a bounding box given constraints in (n-1)
createBoundBox <- function(constr) {
	n <- dim(constr$constr)[2]
	extreme <- findExtremePoints(constr)
	# upper and lower bounds for each dimension in the (n-1) basis
	lb <- apply(extreme, 1, min)
	ub <- apply(extreme, 1, max)
	# bounding function for the length of the hit line
	# x: current position; d: direction of move.
	boundFn <- function(x, d) {
		c(
			max(sapply(1:(n-1), function(i) { min((lb[i] - x[i]) / d[i], (ub[i] - x[i]) / d[i]) })),
			min(sapply(1:(n-1), function(i) { max((lb[i] - x[i]) / d[i], (ub[i] - x[i]) / d[i]) }))
		)
	}
	# starting point (origin)
	start <- (1/(2*(n-1))) * apply(extreme, 1, sum)
	list(bound=boundFn, start=start, lb=lb, ub=ub)
}
