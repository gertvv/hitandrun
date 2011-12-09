# Create a sampler that samples weights from an (n-1) simplex in n-dim space,
# subject to the given constraints.
# n: number of dimensions
# constr: optional matrix of additional constraints
# algo: desired algorithm (default "exactHar")
# startingPoint: NULL or starting point for the chains
# PRECOND: startingPoint is null or has length n
# PRECOND: startingPoint is null or adheres to the constraints
simplex.createSampler <- function(n, constr=NULL, algo="har", startPoint=NULL) {
	stopifnot(is.null(startPoint) || length(startPoint) == n)

	# Initialize (n-1) dimensional sampling space
	basis <- simplex.basis(n)
	a <- simplex.createConstraints(basis, constr)
	hit <- function(x) { TRUE } # we do exact boundary intersections
	bound <- simplex.createBoundBox(a)

	# Set up and verify startPoint
	if (is.null(startPoint)) {
		startPoint <- bound$start
	} else {
		startPoint <- t(basis) %*% (startPoint - (1/n))
	}
	stopifnot(simplex.createHit(a)(startPoint))

	# Create sampling function
	stopifnot(algo %in% c("har", "gibbs"))
	sampler <- NULL
	if (algo == "har") {
		sampler <- function(seed, N) {
			nonAdaptiveHar(seed, N, bound$bound, hit)
		}
	} else {
		sampler <- function(seed, N) {
			gibbs(seed, N, bound$bound, hit)
		}
	}

	# Create chain-generating function
	function(N) {
		result <- sampler(startPoint, N + 1)
		lastPoint <- unlist(tail(result$samples, 1))
		assign("startPoint", lastPoint, envir=parent.env(environment()))
		samples <- simplex.transformResult(basis, tail(result$samples, -1))
		return(list(samples=samples, miss=result$miss))
	}
}

## sample weights from an (n-1) simplex in n-dimensional space, subject to the
## given constraints.
## N: number of desired samples
## n: number of dimensions
## constr: optional matrix of additional constraints
## algo: desired algorithm (default "adaptiveHar")
## startingPoint: NULL or starting point for the chains
## printIters: if > 0, the step size when to print out the current iteration
## PRECOND: startingPoint is null or has length n
## PRECOND: startingPoint is null or adheres to the constraints
simplex.sample <- function(N, n, constr=NULL, algo="adaptiveHar", startingPoint = NULL, printIters = -1) {
  stopifnot(is.null(startingPoint) || length(startingPoint) == n)
  stopifnot(is.null(startingPoint) || all(constr[, 1:n] %*% startingPoint <= constr[, (n + 1)]))
  if (algo == "simplex") { # sample in n-dim space
    hit <- function(x) { TRUE }
    if (!is.null(constr)) {
      hit <- function(x) {
        x <- c(x, -1) # add RHS
        min(constr %*% x <= 0)
      }
    }
    result <- simplexReject(N, n, hit)
    samples <- t(sapply(result$samples, function(x) { x }))
    return(list(samples=samples, miss=result$miss))
  } else { # sample in (n-1)-dim space
    basis <- simplex.basis(n)
    a <- simplex.createConstraints(basis, constr)
    hit <- function(x) { TRUE }
    bound <- simplex.createBoundBox(a)
    if (!is.null(startingPoint)) {
      startingPoint <- t(basis) %*% (startingPoint - (1/n))
    }
    result <- NULL
    if (algo == "adaptiveHar") {
      if (is.null(startingPoint)) {
        result <- adaptiveHar(bound$start, N, bound$bound, hit, printIters)
      }
      else {
        result <- adaptiveHar(startingPoint, N, bound$bound, hit, printIters)
      }
    }
    else if (algo == "nonAdaptiveHar") {
      if (is.null(startingPoint)) {
        result <- nonAdaptiveHar(bound$start, N, bound$bound, hit, printIters)
      }
      else {
        result <- nonAdaptiveHar(startingPoint, N, bound$bound, hit, printIters)
      }
    }
    else if (algo == "gibbs") {
      if (is.null(startingPoint)) {
        result <- gibbs(bound$start, N, bound$bound, hit, printIters)
      }
      else {
        result <- gibbs(startingPoint, N, bound$bound, hit, printIters)
      }
    }
    else if (algo == "bound") {
      result <- boundingBoxReject(N, bound$lb, bound$ub, hit, printIters)
    }
    samples <- simplex.transformResult(basis, result$samples)
    return(list(samples=samples, miss=result$miss))
  }
  stop(paste("Algorithm", algo, "not supported"))
}
