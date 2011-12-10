library(hitandrun)

ordinalConstraints <- function(dim) {
	mergeConstraints(lapply(1:(dim-1), function(i) {ordinalConstraint(dim, i, i + 1)}))
}

n <- 3
transform <- simplex.createTransform(n)
constr <- simplex.createConstraints(transform, ordinalConstraints(n))
seedPoint <- createSeedPoint(constr, homogeneous=TRUE)

N <- 10000
samples <- har(seedPoint, constr, N, 1, homogeneous=TRUE, transform=transform)$samples

# Check dimension
stopifnot(dim(samples) == c(N, n))

# Check that w_i >= w_i+1
stopifnot(sapply(1:(n-1), function(i) {
	all(samples[,i]>=samples[,i+1])
}))

# Check that w_i >= 0
stopifnot(samples >= 0)

# Check that sum_i w_i = 1
E <- 1E-12
stopifnot(apply(samples, 1, sum) > 1 - E)
stopifnot(apply(samples, 1, sum) < 1 + E)

# Check that the points are not all identical
stopifnot(apply(samples, 2, sd) > E)

# Check that seed point is not included in sample
stopifnot(samples[1,] != transform %*% seedPoint)

samples <- har(seedPoint, constr, N, 1, homogeneous=TRUE)$samples

# Check dimension
stopifnot(dim(samples) == c(N, n))

# Check homogeneous coordinate
stopifnot(samples[,n] == 1)

# Check that seed point is not included in sample
stopifnot(samples[1,] != transform %*% seedPoint)
