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
    hit <- simplex.createHit(a)
    bound <- createBoundBox(a)
    if (!is.null(startingPoint)) {
      startingPoint = t(basis) %*% (startingPoint - (1/n))
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
