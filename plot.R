library('rgl')

# sample N sets of n numbers, each number in [0, 1] and each set of numbers summing to 1
# results in a matrix of N rows and n columns
sampleUniform <- function(n, N) {
	sample <- function(n) {
		r <- c(0, sort(runif(n-1)), 1)
		w <- c()
		for (i in 1:n) {
			w[i] <- r[i+1] - r[i]
		}
		w
	}
	t(sapply(1:N, function(x) { sample(n) }))
}

plotWeights <- function(W) {
	plot3d(W, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), xlab="x", ylab="y", zlab="z", axes=F)
	triangles3d(c(0,0,1), c(0,1,0), c(1,0,0), alpha=0.3)
}

## (MOVE TO OTHER FILE)
library('rgl')
library('far') # for orthonormalization
library('lpSolve') # to solve 2(n-1) LPs

drawVector <- function(v, color="black", alpha=0.6) {
	lines3d(c(0, v[1]), c(0, v[2]), c(0, v[3]), color=color, alpha=alpha)
}

testTransformation <- function() {
	# draw the polytope translated to origin
	plot3d(c(0), c(0), c(0), xlim=c(-1/3, 2/3), ylim=c(-1/3, 2/3), zlim=c(-1/3, 2/3))
	triangles3d(c(2/3,-1/3,-1/3), c(-1/3,2/3,-1/3), c(-1/3,-1/3,2/3), alpha=0.3)

	# draw the basis
	basis <- matrix(c(1,0,-1,0,1,-1), ncol=2)
	drawVector((1/3)*basis[,1])
	drawVector((1/3)*basis[,2])

	# draw the orthonormal basis
	on <- orthonormalization(basis, basis=F, norm=T)
	drawVector((1/3)*on[,1], color="red", alpha=1)
	drawVector((1/3)*on[,2], color="red", alpha=1)


	# sample random directions
	dir <- runif(100, 0, 2*pi)
	coord <- matrix(c(cos(dir), sin(dir)), ncol=2)
	samples <- t(sapply(1:dim(coord)[1], function(i) { on %*% coord[i,] }))
	points3d(samples)
}

# HAR

distance2 <- function(x) {
  sum <-  0;
  for (e in x) {
    sum <- sum + (e*e);
  }
  sqrt(sum);
}

normalize2 <- function(x) {
  radius <- distance2(x);
  x/radius;
}  

sampleHyperSphere <- function(dim) {
  x <- rnorm(dim, mean=0, sd=1);
  normalize2(x);
}

# generate a random direction in n-dimensional space
randDir <- function(n) {
	if (n == 2) {
		d <- runif(1, 0, 2*pi)[1] # random direction
		c(cos(d), sin(d)) # transform to x/y
	} else {
		sampleHyperSphere(as.numeric(n))
	}
}

# x0: starting point; niter: number of iterations; bound: bounding box function; hit: hit determination function.
har <- function(x0, niter, bound, hit) {
	stopifnot(as.logical(hit(x0)))
	misses <- 0
	x <- as.list(rep(0, niter))
	n <- length(x0)
	x[[1]] <- x0
	for (i in 2:niter) {
		d <- randDir(n) # random direction
		bounds <- bound(x[[i-1]], d)
		# print(paste("Iteration", i, "x=", paste(x[[i-1]], collapse=","), "d=", paste(d, collapse=",")))
		wasHit <- FALSE
		while (!wasHit) {
			l <- runif(1, bounds[1], bounds[2])[1] # random distance
			xN <- x[[i-1]] + l * d
			# print(paste("b=", paste(bounds, collapse=",")))
			# print(paste("l=", paste(l, collapse=","), "xN=", paste(xN, collapse=",")))
			if (hit(xN)) {
				x[[i]] <- xN
				wasHit <- TRUE
				# print("HIT")
			} else {
				misses <- misses + 1
				x[[i]] <- x[[i-1]]
				# update bounds
				if (l > 0) {
					bounds[2] <- l
				}
				if (l < 0) {
					bounds[1] <- l
				}
				# print(paste("MISS ", i))
			}
		}
	}
	list(x, misses)
}

# translate the n-dimensional constraints to the (n-1)-dimensional space
# basis: orthonormal basis for the (n-1)-dimensional simplex in n-dimensional space
# constr: additional constraints
createConstraintFormulation <- function(basis, userConstr=NULL) {
	n <- dim(basis)[1]

	# basic constraints defining the n-dimensional simplex
	constr1 <- diag(rep(-1, n)) # -1*w[i] <= 0
	rhs1 <- rep(0, n)
	constr2 <- 1 + constr1 # (\sum w[j]) - w[i] <= 1
	rhs2 <- rep(1, n)
	constr <- rbind(constr1, constr2)
	rhs <- c(rhs1, rhs2)

	# user constraints
	if (!is.null(constr)) {
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
createHitSimplex <- function(constr) {
	a <- cbind(a$constr, a$rhs)
	function(x) {
		x <- c(x, 1, -1) # add homogenous coordinate and RHS
		min(a %*% x <= 0)
	}
}

# give the inverse of an invertible matrix
inverse <- function(mat) {
	d <- svd(mat)
	d$v%*%diag(1/d$d)%*%t(d$u)
}

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
	list(bound=boundFn, start=start)
}

# create basis for (translated) n-dim simplex
basis <- function(n) {
	b <- rbind(diag(n-1), rep(-1, n-1))
	orthonormalization(b, basis=F, norm=T)
}

transformResult <- function(basis, samples) {
	n <- dim(basis)[1]
	t(sapply(samples, function(x) { basis %*% x })) + 1/n
}

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

n <- 3
on <- basis(n)
a <- createConstraintFormulation(on)
#a <- createConstraintFormulation(on, rbind(upperRatioConstraint(n, 3, 1, 1.2), lowerRatioConstraint(n, 3, 1, 1/1.2), lowerBoundConstraint(n, 2, 0.2)))
hit <- createHitSimplex(a)
bound <- createBoundBox(a)
samples <- har(bound$start, 1000, bound$bound, hit)
result <- transformResult(on, samples[[1]])
