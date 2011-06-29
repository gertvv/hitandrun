library('rgl')

simplex.plot <- function(W) {
	plot3d(W, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), xlab="x", ylab="y", zlab="z", axes=F)
	triangles3d(c(0,0,1), c(0,1,0), c(1,0,0), alpha=0.3)
}

source('sampleSimplex.R')
source('constraint.R')

n <- 3
constr <- rbind(upperRatioConstraint(n, 3, 1, 1.2), lowerRatioConstraint(n, 3, 1, 1/1.2), lowerBoundConstraint(n, 2, 0.2))
#result <- simplex.sample(1000, n, constr, "simplex")
