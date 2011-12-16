library('rgl')

drawVector <- function(v, color="black", alpha=0.6) {
	lines3d(c(0, v[1]), c(0, v[2]), c(0, v[3]), color=color, alpha=alpha)
}

testTransformation <- function() {
	# draw the polytope translated to origin
	plot3d(c(0), c(0), c(0), xlim=c(-1/3, 2/3), ylim=c(-1/3, 2/3), zlim=c(-1/3, 2/3))
	triangles3d(c(2/3,-1/3,-1/3), c(-1/3,2/3,-1/3), c(-1/3,-1/3,2/3), alpha=0.3)

	# draw the basis
	basis <- matrix(c(1,0,-1,0,1,-1), ncol=2)
	drawVector((1/3)*basis[,1])
	drawVector((1/3)*basis[,2])

	# draw the orthonormal basis
	on <- qr.Q(qr(basis))
	drawVector((1/3)*on[,1], color="red", alpha=1)
	drawVector((1/3)*on[,2], color="red", alpha=1)


	# sample random directions
	dir <- runif(100, 0, 2*pi)
	coord <- matrix(c(cos(dir), sin(dir)), ncol=2)
	samples <- t(sapply(1:dim(coord)[1], function(i) { on %*% coord[i,] }))
	points3d(samples)
}
