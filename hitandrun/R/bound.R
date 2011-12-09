findExtremePoints <- function(constr, homogeneous=FALSE) {
	n <- dim(constr$constr)[2]
	if (homogeneous == TRUE) {
		n <- n - 1
		# add the constraint that the homogenous coordinate equals 1
		constr$constr <- rbind(constr$constr, c(rep(0, n), 1))
		constr$rhs <- c(constr$rhs, 1)
		constr$dir <- c(constr$dir, "=")
	}

	# because lp_solve assumes vars to be non-negative, split each dimension:
	# x_i = y_{2i-1} - y_{2i}
	# thus each dimension is coded as the difference of two variables, give them
	obj <- function(i) {
		obj <- rep(0, 2 * n)
		obj[2 * i - 1] <- 1
		obj[2 * i] <- -1
		obj
	}
	unsplit <- t(sapply(1:n, obj))
	mat <- unsplit
	if (homogeneous == TRUE) { # add row to preserve homogeneous coordinate
		mat <- rbind(unsplit, c(rep(0, 2 * n - 1), 1))
	}
	a <- constr$constr %*% mat

	# for each of the n dimensions, solve 2 LPs to find the min/max
	findExtreme <- function(dir) {
		function(i) {
			lp(dir, obj(i), a, constr$dir, constr$rhs)$solution
		}
	}
	rawExtremes <- cbind(sapply(1:n, findExtreme("min")), sapply(1:n, findExtreme("max")))
	unsplit %*% rawExtremes
}

# generate seed point from constraints (TODO: implement multiple methods)
createSeedPoint <- function(constr, homogeneous=FALSE) {
	n <- dim(constr$constr)[2]
	if (homogeneous == TRUE) {
		n <- n - 1
	}

	extreme <- findExtremePoints(constr, homogeneous)
	# starting point (approximation of centroid)
	p <- (1/(2*n)) * apply(extreme, 1, sum)
	if (homogeneous == TRUE) {
		p[n + 1] <- 1.0 # eliminate floating point imprecision
	}
	p
}

# create a bounding box given constraints in (n-1)
createBoundBox <- function(constr, homogeneous=FALSE) {
	n <- dim(constr$constr)[2]
	extreme <- findExtremePoints(constr, homogeneous)
	# upper and lower bounds for each dimension in the (n-1) basis
	lb <- apply(extreme, 1, min)
	ub <- apply(extreme, 1, max)
	list(lb=lb, ub=ub, extreme=extreme)
}
