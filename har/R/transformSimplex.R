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
		min(a %*% x <= 0)
	}
}
