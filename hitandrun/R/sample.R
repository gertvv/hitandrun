har <- function(x0, constr, N, thin=1, homogeneous=FALSE, transform=NULL) {
	n <- length(x0)
	m <- dim(constr$constr)[1]

	# Verify preconditions that the C-code cannot check
	stopifnot(n == dim(constr$constr)[2])
	stopifnot(m == length(constr$rhs))
	stopifnot(constr$dir == "<=")

	if (homogeneous == FALSE) { # Change to homogeneous coordinates
		n <- n + 1
		constr$constr <- cbind(constr$constr, 0)
		x0 <- c(x0, 1.0)
	}

	samples <- .C("har",
		as.integer(n - 1), as.double(x0),
		as.integer(m), as.double(constr$constr), as.double(constr$rhs),
		as.integer(N), as.integer(1),
		samples=matrix(0.0, nrow=N, ncol=n),
		NAOK=FALSE, DUP=FALSE, PACKAGE="hitandrun"
	)$samples
	if (!is.null(transform)) {
		if (homogeneous == FALSE) { # Add column to eliminate hom. coord.
			transform <- cbind(transform, 0)
		}
		samples <- samples %*% t(transform)
	} else if (homogeneous == FALSE) { # Eliminate hom. coord.
		samples <- samples[,1:(n-1)]
	}
	samples
}

# x0: starting point; niter: number of iterations; bound: bounding box function; hit: hit determination function.
adaptiveHar <- function(x0, niter, bound, hit, printIters=-1) {
	stopifnot(as.logical(hit(x0)))
	misses <- 0
	x <- as.list(rep(0, niter))
	n <- length(x0)
	x[[1]] <- x0
	for (i in 2:niter) {
          if (printIters > 0 && i %% printIters == 0) {
            cat(i, fill = TRUE)
          }

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

# x0: starting point; niter: number of iterations; bound: bounding box function; hit: hit determination function.
nonAdaptiveHar <- function(x0, niter, bound, hit, printIters=-1) {
	stopifnot(as.logical(hit(x0)))
	misses <- 0
	x <- as.list(rep(0, niter))
	n <- length(x0)
	x[[1]] <- x0
	for (i in 2:niter) {
          if (printIters > 0 && i %% printIters == 0) {
            cat(i, fill = TRUE)
          }

		d <- randDir(n) # random direction
		bounds <- bound(x[[i-1]], d)
		l <- runif(1, bounds[1], bounds[2])[1] # random distance
		xN <- x[[i-1]] + l * d
		if (hit(xN)) {
			x[[i]] <- xN
		} else {
			x[[i]] <- x[[i-1]]
			misses <- misses + 1
		}
	}
	list(samples=x, miss=misses)
}

# x0: starting point; niter: number of iterations; bound: bounding box function; hit: hit determination function.
gibbs <- function(x0, niter, bound, hit, printIters=-1) {
	stopifnot(as.logical(hit(x0)))
	misses <- 0
	x <- as.list(rep(0, niter))
	n <- length(x0)
	x[[1]] <- x0
	# print(x0)
	for (i in 2:niter) {
          if (printIters > 0 && i %% printIters == 0) {
            cat(i, fill = TRUE)
          }          
		xN <- x[[i - 1]]
		for (j in 1:n) {
			d <- rep(0, n)
			d[j] <- 1 # update the i-th dimension
			# print(paste("x=", xN))
			# print(paste("d=", d))
			bounds <- bound(xN, d)
			# print(paste("b=", paste(bounds, collapse=",")))
			l <- runif(1, bounds[1], bounds[2])[1] # random distance
			xM <- xN + l * d
			if (hit(xM)) {
				xN <- xM
			} else {
				misses <- misses + 1
			}
		}
		x[[i]] <- xN
	}
	list(samples=x, miss=misses)
}

# bounding box rejection sampling
boundingBoxReject <- function(niter, lb, ub, hit, printIters=-1) {
	misses <- 0
	x <- as.list(rep(0, niter))
	n <- length(lb)
	d <- ub - lb # pre-calculate for use in loop
	for (i in 1:niter) {
          if (printIters > 0 && i %% printIters == 0) {
            cat(i, fill = TRUE)
          }          
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
simplexReject <- function(niter, n, hit, printIters=-1) {
	misses <- 0
	x <- as.list(rep(0, niter))
	for (i in 1:niter) {
          if (printIters > 0 && i %% printIters == 0) {
            cat(i, fill = TRUE)
          }          
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
