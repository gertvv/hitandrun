hitandrun <- function(constr,
    n.samples=1E4,
    thin.fn = function(n) { ceiling(log(n + 1)/4 * n^3) },
    thin = NULL,
    x0 = NULL) {
  stopifnot(length(constr[['rhs']]) == length(constr[['dir']]))
  stopifnot(length(constr[['rhs']]) == nrow(constr[['constr']]))
  stopifnot(length(constr[['rhs']]) > 0)

  # Separate equality from inequality constraints
  sel.eq <- constr$dir == '='
  eq.constr <- filterConstraints(constr, sel.eq)
  ineq.constr <- filterConstraints(constr, !sel.eq)
  stopifnot(all(ineq.constr[['dir']] == '<='))

  # Generate basis, transform, seed point
  basis <- 
  if (length(eq.constr$dir) > 0) {
    solution.basis(eq.constr)
  } else {
    n <- ncol(ineq.constr$constr)
    list(translate=rep(0, n), basis=diag(n))
  }
  
  transform <- createTransform(basis)
  constr <- transformConstraints(transform, ineq.constr)
  if (is.null(x0)) {
    x0 <- createSeedPoint(constr, homogeneous=TRUE)
  }

  # sample
  n <- length(x0) - 1
  if (n == 0) {
    list(samples = matrix(rep(basis$translate, each=n.samples), nrow=n.samples),
         xN = 1)
  } else {
    if (is.null(thin)) {
      thin <- thin.fn(n)
    }
    har(x0, constr, N=n.samples * thin, thin=thin, homogeneous=TRUE, transform=transform)
  }
}
