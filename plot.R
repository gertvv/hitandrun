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

drawVector <- function(v, color="black", alpha=0.6) {
	lines3d(c(0, v[1]), c(0, v[2]), c(0, v[3]), color=color, alpha=alpha)
}

testTransformation <- function() {
	# draw the polytope translated to origin
	plot3d(c(0), c(0), c(0), xlim=c(-1/3, 2/3), ylim=c(-1/3, 2/3), zlim=c(-1/3, 2/3))
	triangles3d(c(2/3,-1/3,-1/3), c(-1/3,2/3,-1/3), c(-1/3,-1/3,2/3), alpha=0.3)

	# draw the basis
	drawVector((1/3)*basis[,1])
	drawVector((1/3)*basis[,2])

	# draw the orthonormal basis
	drawVector((1/3)*on[,1], color="red", alpha=1)
	drawVector((1/3)*on[,2], color="red", alpha=1)


	# sample random directions
	dir <- runif(100, 0, 2*pi)
	coord <- matrix(c(cos(dir), sin(dir)), ncol=2)
	samples <- t(sapply(1:dim(coord)[1], function(i) { on %*% coord[i,] }))
	points3d(samples)
}

# HAR

# generate a random direction in n-dimensional space
# FIXME: implement for n != 2
randDir <- function(n) {
	stopifnot(n == 2)
	d <- runif(1, 0, 2*pi)[1] # random direction
	c(cos(d), sin(d)) # transform to x/y
}

# x0: starting point; niter: number of iterations; bound: bounding box function; hit: hit determination function.
har <- function(x0, niter, bound, hit) {
	misses <- 0
	x <- as.list(rep(0, niter))
	n <- dim(x0)
	x[[1]] <- x0
	for (i in 2:niter) {
		d <- randDir(n) # random direction
		bounds <- bound(x[[i-1]], d)
		wasHit <- FALSE
		while (!wasHit) {
			l <- runif(1, bounds[1], bounds[2])[1] # random distance
			xN <- x[[i-1]] + l * d
			if (hit(xN)) {
				x[[i]] <- xN
				wasHit <- TRUE
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
			}
		}
	}
	list(x, misses)
}

# hit calculation for (n-1)-dimensional simplex
# basis: orthonormal basis for the (n-1)-dimensional simplex in n-dimensional space
createHitSimplex <- function(basis) {
	n <- dim(basis)[1]
	a1 <- diag(rep(-1, n)) # -1*w[i] <= 0
	a2 <- 1 + c1 # (\sum w[j]) - w[i] <= 1
	a <- cbind(rbind(a1, a2), c(rep(0, n), rep(1, n)))
	function(x) {
		x <- basis %*% x # change of basis
		w <- rbind(1/n + x, -1) # translation
		min(a %*% w <= 0)
	}
}

# give the inverse of an invertible matrix
inverse <- function(mat) {
	d <- svd(mat)
	d$v%*%diag(1/d$d)%*%t(d$u)
}

# create a bounding box given a (n-1) basis and extreme points in R^n
createBoundBox <- function(basis, extreme) {
	n <- dim(basis)[2]
	extreme <- t(basis) %*% (extreme - 1/n) # extreme points
	lb <- apply(extreme, 1, min)
	ub <- apply(extreme, 1, max)
	function(x, d) {
		c(
			max(sapply(1:n, function(i) { min((lb[i] - x[i]) / d[i], (ub[i] - x[i]) / d[i]) })),
			min(sapply(1:n, function(i) { max((lb[i] - x[i]) / d[i], (ub[i] - x[i]) / d[i]) }))
		)
	}
}

basis <- matrix(c(1, 0, -1, 0, 1, -1), nrow=3)
on <- orthonormalization(basis, basis=F, norm=T)
hit <- createHitSimplex(on)
bound <- createBoundBox(on, diag(3))
samples <- har(c(0,0), 1000, bound, hit)
result <- t(sapply(samples[[1]], function(x) { on %*% x })) + 1/3

