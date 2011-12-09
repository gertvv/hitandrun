# create basis for (translated) n-dim simplex
simplex.basis <- function(n) {
	b <- rbind(diag(n-1), rep(-1, n-1))
	orthonormalization(b, basis=F, norm=T)
}

# Generate a projection matrix that transforms an (n-1) dimensional vector in
# homogeneous coordinate representation to an n-dimensional weight vector.
simplex.createTransform <- function(n) {
	basis <- simplex.basis(n)
	# add one extra element to vectors in each basis (homogeneous coordinate
	# representation)
	basis <- rbind(cbind(basis, rep(0, n)), c(rep(0, n - 1), 1))
	# create translation matrix (using homogenous coordinates)
	translation <- cbind(diag(n), rep(1/n, n))
	# successively apply basis transformation and translation
	translation %*% basis
}

# translate the n-dimensional constraints to the (n-1)-dimensional space
# transform: transform created by simplex.createTransform 
# userConstr: additional constraints
simplex.createConstraints <- function(transform, userConstr=NULL) {
	n <- dim(transform)[1]

	# basic constraints defining the (n-1)-dimensional simplex
	constr <- diag(rep(-1, n)) # -1*w[i] <= 0
	rhs <- rep(0, n)

	# user constraints
	if (!is.null(userConstr)) {
		stopifnot(userConstr$dir == "<=")
		constr <- rbind(constr, userConstr$constr)
		rhs <- c(rhs, userConstr$rhs)
	}

	constr <- constr %*% transform

	# give directions
	m <- dim(constr)[1]
	dir <- rep("<=", m)

	list(constr=constr, rhs=rhs, dir=dir)
}
