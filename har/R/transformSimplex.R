# create basis for (translated) n-dim simplex
simplex.basis <- function(n) {
	b <- rbind(diag(n-1), rep(-1, n-1))
	orthonormalization(b, basis=F, norm=T)
}

# transform results back from (n-1)-dim simplex space to n-dim space
simplex.transformResult <- function(basis, samples) {
	n <- dim(basis)[1]
	t(sapply(samples, function(x) { basis %*% x })) + 1/n
}

# translate the n-dimensional constraints to the (n-1)-dimensional space
# basis: orthonormal basis for the (n-1)-dimensional simplex in n-dimensional space
# constr: additional constraints
simplex.createConstraints <- function(basis, userConstr=NULL) {
	n <- dim(basis)[1]

	# basic constraints defining the (n-1)-dimensional simplex
	constr <- diag(rep(-1, n)) # -1*w[i] <= 0
	rhs <- rep(0, n)

	# user constraints
	if (!is.null(userConstr)) {
		constr <- rbind(constr, userConstr[, 1:n])
		rhs <- c(rhs, userConstr[, n + 1])
	}

	# add one extra element to vectors in each basis (homogenous coordinate representation)
	basis <- rbind(cbind(basis, rep(0, n)), c(rep(0, n - 1), 1))
	# create translation matrix (using homogenous coordinates)
	translation <- rbind(cbind(diag(n), rep(1/n, n)), c(rep(0, n), 1))
	# successively apply basis transformation, translation and then constraint calculation
	constr <- cbind(constr, 0) %*% translation %*% basis

	# add the constraint that the homogenous coordinate equals 1
	m <- dim(constr)[1]
	constr <- rbind(constr, c(rep(0, n-1), 1))
	rhs <- c(rhs, 1)

	# give directions
	dir <- c(rep("<=", m), "=")

	list(constr=constr, rhs=rhs, dir=dir)
}

# hit calculation for (n-1)-dimensional simplex
# constr: constraint formulation in (n-1)-dimensional space
simplex.createHit <- function(constr) {
	a <- cbind(constr$constr, constr$rhs)
	function(x) {
		x <- c(x, 1, -1) # add homogenous coordinate and RHS
		all(a %*% x <= 0)
	}
}

# create a bounding box given constraints in (n-1)
simplex.createBoundBox <- function(constr) {
	n <- dim(constr$constr)[2]
	extreme <- findExtremePoints(constr)
	# upper and lower bounds for each dimension in the (n-1) basis
	lb <- apply(extreme, 1, min)
	ub <- apply(extreme, 1, max)
	# bounding function for the length of the hit line
	# x: current position; d: direction of move.
	boundFn <- function(x, d) {
		x <- c(x, 1) # add homogeneous coordinate
		d <- c(d, 0) # add homogeneous coordinate
		# we know Ax <= b, now we need to find the values of t such that
		# A(x + td) <= b, i.e. t(Ad) <= b - Ax:
		# T = [
		#   max_{i:(Ad)_i<0} (b - Ax)_i / (Ad)_i,
		#   min_{i:(Ad)_i>0} (b - Ax)_i / (Ad)_i
		# ]
		a <- constr$rhs - constr$constr %*% x #(b - Ax)
		c <- constr$constr %*% d #Ad
		c(
			max(a[c < 0] / c[c < 0]),
			min(a[c > 0] / c[c > 0])
		)
	}
	# starting point (approximation of centroid)
	start <- (1/(2*(n-1))) * apply(extreme, 1, sum)
	list(bound=boundFn, start=start, lb=lb, ub=ub)
}
