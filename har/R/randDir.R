distance2 <- function(x) {
  sum <-  0;
  for (e in x) {
    sum <- sum + (e*e);
  }
  sqrt(sum);
}

normalize2 <- function(x) {
  radius <- distance2(x);
  x/radius;
}  

sampleHyperSphere <- function(dim) {
  x <- rnorm(dim, mean=0, sd=1);
  normalize2(x);
}

# generate a random direction in n-dimensional space
randDir <- function(n) {
	if (n == 2) {
		d <- runif(1, 0, 2*pi)[1] # random direction
		c(cos(d), sin(d)) # transform to x/y
	} else {
		sampleHyperSphere(as.numeric(n))
	}
}
