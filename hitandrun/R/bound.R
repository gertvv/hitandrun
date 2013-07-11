homogeneousCoordinateConstraint <- function(n) {
	list(constr=c(rep(0, n), 1), rhs=c(1), dir=c("="))
}

findExtremePoints <- function(constr, homogeneous=FALSE) {
	n <- dim(constr$constr)[2]
	nh <- n
	h <- if (homogeneous == TRUE) {
		n <- n - 1
		hom <- homogeneousCoordinateConstraint(n)
		makeH(constr$constr, constr$rhs, hom$constr, hom$rhs)
	} else {
		makeH(constr$constr, constr$rhs)
	}

	# for each of the n dimensions, solve 2 LPs to find the min/max
	findExtreme <- function(minimize) {
		function(i) {
			obj <- rep(0, nh)
			obj[i] <- 1
			lpcdd(h, obj, minimize=minimize)$primal.solution
		}
	}
	t(cbind(sapply(1:n, findExtreme(TRUE)), sapply(1:n, findExtreme(FALSE))))
}

findVertices <- function(constr, homogeneous=FALSE) {
	h <- if (homogeneous == TRUE) {
		n <- dim(constr$constr)[2]
		hom <- homogeneousCoordinateConstraint(n - 1)
		makeH(constr$constr, constr$rhs, hom$constr, hom$rhs)
	} else {
		makeH(constr$constr, constr$rhs)
	}

	v <- q2d(scdd(d2q(h))$output)
	# Check that the output are vertices, not other things that would indicate
	# a bug
	stopifnot(v[,1] == "0")
	stopifnot(v[,2] == "1")

	# Return the vertices only 
	v[,-c(1,2)]
}

# Finds points on the boundary of the polytope so that the returned set of points spans the polytope.
# Assumes the polytope is of equal dimension to the space.
findBoundaryPoints <- function(constr, homogeneous=FALSE) {
	n <- ncol(constr$constr) - homogeneous # dimensionality of space
	extreme <- findExtremePoints(constr, homogeneous)
	m <- qr(extreme[,1:n])$rank
	while (m < n) {
		B <- qr.Q(qr(t(extreme[,1:n])))
		if (homogeneous) {
			B <- rbind(cbind(B, rep(0, n)), c(rep(0, n), 1))
		}
		extreme <- findExtremePoints(
			list(constr = constr$constr %*% B, dir = constr$dir, rhs = constr$rhs),
			homogeneous = homogeneous) %*% t(B)
		m1 <- qr(extreme[,1:n])$rank
		stopifnot(m1 > m)
		m <- m1
	}
	extreme
}

# generate seed point from constraints
createSeedPoint <- function(constr, homogeneous=FALSE, randomize=FALSE,
		method="extremes") {
	stopifnot(method %in% c("extremes", "vertices"))

	n <- dim(constr$constr)[2]
	if (homogeneous == TRUE) {
		n <- n - 1
	}

	extreme <- if (method == "extremes") {
		findBoundaryPoints(constr, homogeneous)
	} else {
		findVertices(constr, homogeneous)
	}

	# starting point 
	m <- dim(extreme)[1]
	p <- if (randomize == TRUE) { # random weighting
		w <- as.vector(simplex.sample(m, 1)$samples)
		apply(extreme, 2, function(row) { sum(w * row) })
	} else { # mean: approximation of centroid
		(1/m) * apply(extreme, 2, sum)
	}

	if (homogeneous == TRUE) {
		p[n + 1] <- 1.0 # eliminate floating point imprecision
	}
	p
}

# create a bounding box given constraints in (n-1)
createBoundBox <- function(constr, homogeneous=FALSE) {
	n <- dim(constr$constr)[2]
	extreme <- findExtremePoints(constr, homogeneous)
	if (homogeneous == TRUE) {
		extreme <- extreme[ , 1:(n-1)]
	}
	# upper and lower bounds for each dimension in the (n-1) basis
	lb <- apply(extreme, 2, min)
	ub <- apply(extreme, 2, max)
	list(lb=lb, ub=ub)
}
