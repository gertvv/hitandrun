eq.constr <- function(constr) {
  filterConstraints(constr, constr[['dir']] == '=')
}

iq.constr <- function(constr) {
  x <- filterConstraints(constr, constr[['dir']] != '=')
  stopifnot(all(x[['dir']] == '<='))
  x
}

eliminateRedundant <- function(constr) {
  n <- ncol(constr$constr)

  eq <- eq.constr(constr)
  iq <- iq.constr(constr)

  h <- 
    if (length(eq[['dir']]) > 0) {
      rcdd::makeH(iq$constr, iq$rhs, eq$constr, eq$rhs)
    } else {
      rcdd::makeH(iq$constr, iq$rhs)
    }

  h.nr <- rcdd::q2d(rcdd::redundant(rcdd::d2q(h))$output)
  rows <- h.nr[, 1] == 0
  nvar <- ncol(constr$constr)
  list(
    constr = -h.nr[ , -c(1, 2), drop=FALSE],
    rhs    = h.nr[ , 2, drop=TRUE],
    dir    = c("<=", "=")[h.nr[, 1] + 1])
}

har.init <- function(constr,
    thin.fn = function(n) { ceiling(log(n + 1)/4 * n^3) },
    thin = NULL,
    x0.randomize = FALSE, x0.method="slacklp",
    x0 = NULL) {
  stopifnot(length(constr[['rhs']]) == length(constr[['dir']]))
  stopifnot(length(constr[['rhs']]) == nrow(constr[['constr']]))
  stopifnot(length(constr[['rhs']]) > 0)

  # Eliminate redundant constraints to catch any implied equalities
  constr <- eliminateRedundant(constr)

  # Separate equality from inequality constraints
  eq <- eq.constr(constr)
  iq <- iq.constr(constr)

  # Generate basis, transform, seed point
  basis <- 
  if (length(eq$dir) > 0) {
    solution.basis(eq)
  } else {
    n <- ncol(constr$constr)
    list(translate=rep(0, n), basis=diag(n))
  }
  
  transform <- createTransform(basis)
  constr <- transformConstraints(transform, iq)
  if (is.null(x0)) {
    x0 <- createSeedPoint(constr, homogeneous=TRUE, randomize=x0.randomize, method=x0.method)
  }

  n <- length(x0) - 1
  if (is.null(thin)) {
    thin <- if (n == 0) 1 else thin.fn(n)
  }

  list(
    basis = basis,
    transform = transform,
    constr = constr,
    x0 = x0,
    thin = thin)
}

har.run <- function(state, n.samples) {
  result <- with(state, {
    n <- length(x0) - 1
    if (n == 0) {
      list(samples = matrix(rep(basis$translate, each=n.samples), nrow=n.samples), xN = 1)
    } else {
      har(x0, constr, N=n.samples * thin, thin=thin, homogeneous=TRUE, transform=transform)
    }
  })
  state$x0 <- result$xN
  list(state = state, samples = result$samples)
}

hitandrun <- function(constr,
    n.samples=1E4,
    thin.fn = function(n) { ceiling(log(n + 1)/4 * n^3) },
    thin = NULL,
    x0.randomize = FALSE, x0.method="slacklp",
    x0 = NULL) {
  state <- har.init(constr, thin.fn, thin, x0.randomize, x0.method, x0)
  result <- har.run(state, n.samples)
  result$samples
}
