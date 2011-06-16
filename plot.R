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

# draw the polytope translated to origin
plot3d(c(0), c(0), c(0), xlim=c(-1/3, 2/3), ylim=c(-1/3, 2/3), zlim=c(-1/3, 2/3))
triangles3d(c(2/3,-1/3,-1/3), c(-1/3,2/3,-1/3), c(-1/3,-1/3,2/3), alpha=0.3)

# draw the basis
basis <- matrix(c(1, 0, -1, 0, 1, -1), nrow=3)
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


# HAR

# x0: starting point; niter: number of iterations; bound: bounding box function; hit: hit determination function.
har <- function(x0, niter, bound, hit) {
	x <- as.list(rep(0, niter))
	x[[1]] <- x0
	wasHit <- TRUE
	d <- c(0, 0)
	bounds <- c(0, 0)
	for (i in 2:niter) {
		print(wasHit)
		if (wasHit) {
			d <- runif(1, 0, 2*pi)[1] # random direction
			d <- c(cos(d), sin(d))
			bounds <- bound(x[[i-1]], d)
		}
		print(d)
		print(bounds)
		l <- runif(1, bounds[1], bounds[2])[1] # random distance
		xN <- x[[i]] + l * d
		print(xN)
		if (hit(xN)) {
			x[[i]] <- xN
			wasHit <- TRUE
		} else {
			x[[i]] <- x[[i-1]]
			wasHit <- FALSE
			# update bounds
			if (l > 0) {
				bounds[2] <- l
			}
			if (l < 0) 
				bounds[1] <- l
			}
		}
	}
	x
}

# x0: starting point; niter: number of iterations; bound: bounding box function; hit: hit determination function.
har <- function(x0, niter, bound, hit) {
	x <- as.list(rep(0, niter))
	x[[1]] <- x0
	for (i in 2:niter) {
		print(paste(c("Iteration", i)))
		d <- runif(1, 0, 2*pi)[1] # random direction
		d <- c(cos(d), sin(d)) # transform to x/y
		bounds <- bound(x[[i-1]], d)
		print(paste(c("x=", x[[i-1]], "d=", d, "b=", bounds)))
		wasHit <- FALSE
		while (!wasHit) {
			l <- runif(1, bounds[1], bounds[2])[1] # random distance
			xN <- x[[i-1]] + l * d
			print(paste(c("l=", l, "xN=", xN)))
			if (hit(xN)) {
				x[[i]] <- xN
				wasHit <- TRUE
				print("HIT")
			} else {
				x[[i]] <- x[[i-1]]
				# update bounds
				if (l > 0) {
					bounds[2] <- l
				}
				if (l < 0) {
					bounds[1] <- l
				}
				print(paste(c("MISS; b=", bounds)))
			}
		}
	}
	x
}

# a bounding box of 2/3 around origin
boundBox <- function(x, d) {
	c(
		max(sapply(1:length(x), function(i) { min((-2/3 - x[i]) / d[i], (2/3 - x[i]) / d[i]) })),
		min(sapply(1:length(x), function(i) { max((-2/3 - x[i]) / d[i], (2/3 - x[i]) / d[i]) }))
	)
}

EPSILON <- 0.000001

# hit calculation for (n-1)-dimensional simplex
# basis: orthonormal basis for the (n-1)-dimensional simplex in n-dimensional space
createHitSimplex <- function(basis) {
#	correction <- basis %*% t(basis)
	function(x) {
		x <- basis %*% x # change of basis
#		x <- correction %*% x # project onto plane
		w <- 1/3 + x # translation (FIXME: correct for all n?)
		allBounds <- min(sapply(w, function(wi) { wi >= 0 }))
		sumBounds <- min(c(w[1] + w[2] <= 1, w[1] + w[3] <= 1, w[2] + w[3] <= 1))
		allBounds && sumBounds
	}
}


hit <- createHitSimplex(on)
samples <- har(c(0,0), 1000, boundBox, hit)
result <- t(sapply(samples, function(x) { on %*% x })) + 1/3

