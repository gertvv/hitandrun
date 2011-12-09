# Temporary file to test native implementation of HAR.

library(har)

source('har/tests/simplex3d.R')

library('rgl')

simplex.plot <- function(W) {
        plot3d(W, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), xlab="x", ylab="y", zlab="z", axes=F)
        triangles3d(c(0,0,1), c(0,1,0), c(1,0,0), alpha=0.3)
}

if (n == 3) {
	simplex.plot(samples)
}
