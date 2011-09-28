# sample weights from an (n-1) simplex in n-dimensional space, subject to the
# given constraints.
# N: number of desired samples
# n: number of dimensions
# constr: optional matrix of additional constraints
# algo: desired algorithm (default "adaptiveHar")
simplex.sample <- function(N, n, constr=NULL, algo="adaptiveHar") {
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
		hit <- simplex.createHit(a)
		bound <- createBoundBox(a)
		result <- NULL
		if (algo == "adaptiveHar") {
			result <- adaptiveHar(bound$start, N, bound$bound, hit)
		} else if (algo == "nonAdaptiveHar") {
			result <- nonAdaptiveHar(bound$start, N, bound$bound, hit)
		} else if (algo == "gibbs") {
			result <- gibbs(bound$start, N, bound$bound, hit)
		} else if (algo == "bound") {
			result <- boundingBoxReject(N, bound$lb, bound$ub, hit)
		}
		samples <- simplex.transformResult(basis, result$samples)
		return(list(samples=samples, miss=result$miss))
	}
	stop(paste("Algorithm", algo, "not supported"))
}
