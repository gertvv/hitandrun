# Temporary file to test native implementation of HAR.

library(har)

ordinalConstraints <- function(dim) {
	t(sapply(1:(dim-1), function(i) {ordinalConstraint(dim, i, i + 1)}))
}

library('rgl')

simplex.plot <- function(W) {
        plot3d(W, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), xlab="x", ylab="y", zlab="z", axes=F)
        triangles3d(c(0,0,1), c(0,1,0), c(1,0,0), alpha=0.3)
}

n <- 3
constr <- simplex.createConstraints(basis=simplex.basis(n), userConstr=ordinalConstraints(n))
bb <- simplex.createBoundBox(constr)
m <- length(constr$rhs) - 1

dyn.load("har.so")
N <- 10000
samples <- .C("har", as.integer(n - 1), as.double(c(bb$start, 1.0)), as.integer(m), as.double(constr$constr[1:m,]), as.double(constr$rhs[1:m]), as.integer(N), as.integer(1), samples=matrix(0.0, nrow=N, ncol=n))$samples

basis <- rbind(cbind(simplex.basis(n), rep(0, n)), c(rep(0, n - 1), 1))
translation <- cbind(diag(n), rep(1/n, n))
transform <- translation %*% basis

trSamples <- t(transform %*% t(samples))

# Check that w_i >= w_i+1
stopifnot(sapply(1:(n-1), function(i) {
	all(trSamples[,i]>=trSamples[,i+1])
}))

# Check that w_i >= 0
stopifnot(trSamples >= 0)

# Check that sum_i w_i = 1
E <- 1E-12
stopifnot(apply(trSamples, 1, sum) > 1 - E)
stopifnot(apply(trSamples, 1, sum) < 1 + E)

# Check that the points are not all identical
stopifnot(apply(trSamples, 2, sd) > E)

if (n == 3) {
	simplex.plot(trSamples)
}
