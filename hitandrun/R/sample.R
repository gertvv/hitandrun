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

# bounding box rejection sampling
# TODO: implement in C
bbReject <- function(lb, ub, constr, N, homogeneous=FALSE, transform=NULL) {
	error("MUST FIX")
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
# TODO: implement in C
simplex.sample <- function(n, N) {
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
