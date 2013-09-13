har <- function(x0, constr, N, thin=1, homogeneous=FALSE, transform=NULL) {
  n <- length(x0)
  m <- nrow(constr$constr)

  # Verify preconditions that the C-code cannot check
  stopifnot(n > homogeneous)
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
    as.integer(N), as.integer(thin),
    samples=matrix(0.0, nrow=N/thin, ncol=n),
    NAOK=FALSE, DUP=FALSE, PACKAGE="hitandrun"
  )$samples
  xN <- samples[N/thin, , drop=TRUE]
  if (!is.null(transform)) {
    if (homogeneous == FALSE) { # Add column to eliminate hom. coord.
      transform <- cbind(transform, 0)
    }
    samples <- samples %*% t(transform)
  } else if (homogeneous == FALSE) { # Eliminate hom. coord.
    samples <- samples[ , 1:(n-1), drop=FALSE]
  }
  list(samples=samples, xN=xN)
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

  result <- .C("bbReject",
    as.integer(n - 1), as.double(lb), as.double(ub),
    as.integer(m), as.double(constr$constr), as.double(constr$rhs),
    as.integer(N),
    samples=matrix(0.0, nrow=N, ncol=n), reject=double(1),
    NAOK=FALSE, DUP=FALSE, PACKAGE="hitandrun"
  )
  samples <- result$samples
  if (!is.null(transform)) {
    if (homogeneous == FALSE) { # Add column to eliminate hom. coord.
      transform <- cbind(transform, 0)
    }
    samples <- samples %*% t(transform)
  } else if (homogeneous == FALSE) { # Eliminate hom. coord.
    samples <- samples[, 1:(n-1), drop=FALSE]
  }
  list(samples=samples, rejectionRate=result$reject)
}

simplex.sample <- function(n, N, sort=FALSE) {
  samples <- .C("simplexSample",
    as.integer(n), as.integer(sort), as.integer(N),
    samples=matrix(0.0, nrow=N, ncol=n),
    NAOK=FALSE, DUP=FALSE, PACKAGE="hitandrun"
  )$samples
  list(samples=samples)
}

hypersphere.sample <- function(n, N) { 
  sample <- .C("randDirForR",
               sample=matrix(0.0, nrow=N, ncol=n),
               as.integer(n),
               as.integer(N),
               NAOK=FALSE, DUP=FALSE, PACKAGE="hitandrun"
               )$sample
  
  return(sample)
}
