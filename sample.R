source('randDir.R')

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
	list(samples=x, miss=misses)
}

# bounding box rejection sampling
boundingBoxReject <- function(niter, lb, ub, hit) {
	misses <- 0
	x <- as.list(rep(0, niter))
	n <- length(lb)
	d <- ub - lb # pre-calculate for use in loop
	for (i in 1:niter) {
		wasHit <- FALSE
		while (!wasHit) {
			xN <- lb + d * runif(n)
			if (hit(xN)) {
				x[[i]] <- xN
				wasHit <- TRUE
			} else {
				misses <- misses + 1
			}
		}
	}
	list(samples=x, miss=misses)
}

# sample N sets of n numbers, each number in [0, 1] and each set of numbers summing to 1 (i.e. numbers from the n-simplex)
# results in a matrix of N rows and n columns
simplex.baseSample <- function(n, N) {
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

# simplex rejection sampling
simplexReject <- function(niter, n, hit) {
	misses <- 0
	x <- as.list(rep(0, niter))
	for (i in 1:niter) {
		wasHit <- FALSE
		while (!wasHit) {
			xN <- simplex.baseSample(n, 1)
			if (hit(xN)) {
				x[[i]] <- xN
				wasHit <- TRUE
			} else {
				misses <- misses + 1
			}
		}
	}
	list(samples=x, miss=misses)
}
