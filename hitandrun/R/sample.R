findFace <- function(x, constr) {
  stopifnot(length(x) == ncol(constr$constr))
  d <- constr$constr %*% x - constr$rhs
  which.max(d)
}

har <- function(x0, constr, N, thin=1, homogeneous=FALSE, transform=NULL) {
  n <- length(x0)
  m <- nrow(constr$constr)

  # Verify preconditions
  stopifnot(n > homogeneous)
  stopifnot(n == ncol(constr$constr))
  stopifnot(m == length(constr$rhs))
  stopifnot(constr$dir == "<=")

  if (homogeneous == FALSE) { # Change to homogeneous coordinates
    n <- n + 1
    constr$constr <- cbind(constr$constr, 0)
    x0 <- c(x0, 1.0)
  }

  stopifnot(x0[n] == 1.0)
  stopifnot(N %% thin == 0)

  rval <- .Call("hitandrun_har", x0, constr$constr, constr$rhs, N, thin)
  result <- list(samples=rval)

  result$xN <- result$samples[N/thin, , drop=TRUE]
  if (!is.null(transform)) {
    if (homogeneous == FALSE) { # Add column to eliminate hom. coord.
      transform <- cbind(transform, 0)
    }
    result$samples <- result$samples %*% t(transform)
  } else if (homogeneous == FALSE) { # Eliminate hom. coord.
    result$samples <- result$samples[ , 1:(n-1), drop=FALSE]
  }
  result
}

sab <- function(x0, i0, constr, N, thin=1, homogeneous=FALSE, transform=NULL) {
  n <- length(x0)
  m <- nrow(constr$constr)

  # Verify preconditions
  stopifnot(n > homogeneous)
  stopifnot(n == ncol(constr$constr))
  stopifnot(m == length(constr$rhs))
  stopifnot(constr$dir == "<=")

  if (homogeneous == FALSE) { # Change to homogeneous coordinates
    n <- n + 1
    constr$constr <- cbind(constr$constr, 0)
    x0 <- c(x0, 1.0)
  }

  stopifnot(x0[n] == 1.0)
  stopifnot(N %% thin == 0)

  # normalize the constraints (required for shake-and-bake)
  for (i in 1:m) {
    norm <- sqrt(sum(constr$constr[i,1:(n-1)]^2))
    constr$constr[i,] <- constr$constr[i,] / norm
    constr$rhs[i] <- constr$rhs[i] / norm
  }

  rval <- .Call("hitandrun_sab", x0, i0, constr$constr, constr$rhs, N, thin)
  result <- list(samples=rval[[1]], faces=rval[[2]])

  result$xN <- result$samples[N/thin, , drop=TRUE]
  result$iN <- result$faces[N/thin]
  if (!is.null(transform)) {
    if (homogeneous == FALSE) { # Add column to eliminate hom. coord.
      transform <- cbind(transform, 0)
    }
    result$samples <- result$samples %*% t(transform)
  } else if (homogeneous == FALSE) { # Eliminate hom. coord.
    result$samples <- result$samples[ , 1:(n-1), drop=FALSE]
  }
  result
}

bbReject <- function(lb, ub, constr, N, homogeneous=FALSE, transform=NULL) {
  n <- ncol(constr$constr)
  m <- nrow(constr$constr)

  # Verify preconditions that the C-code cannot check
  stopifnot(n > homogeneous)
  if (homogeneous == FALSE) {
    stopifnot(n == length(lb))
    stopifnot(n == length(ub))
  } else {
    stopifnot(n - 1 == length(lb))
    stopifnot(n - 1 == length(ub))
  }
  stopifnot(m == length(constr$rhs))
  stopifnot(constr$dir == "<=")

  if (homogeneous == FALSE) { # Change to homogeneous coordinates
    n <- n + 1
    constr$constr <- cbind(constr$constr, 0)
  }

  result <- .Call("hitandrun_bbReject", lb, ub, constr$constr, constr$rhs, N)
  samples <- result[[1]]
  reject <- result[[2]]
  if (!is.null(transform)) {
    if (homogeneous == FALSE) { # Add column to eliminate hom. coord.
      transform <- cbind(transform, 0)
    }
    samples <- samples %*% t(transform)
  } else if (homogeneous == FALSE) { # Eliminate hom. coord.
    samples <- samples[, 1:(n-1), drop=FALSE]
  }
  list(samples=samples, rejectionRate=reject)
}

simplex.sample <- function(n, N, sort=FALSE) {
  samples <- .Call("hitandrun_simplexSample", n, sort, N);
  list(samples=samples)
}

hypersphere.sample <- function(n, N) { 
  .Call("hitandrun_hypersphereSample", n, N)
}
