distance2 <- function(x) {
  sum <-  0;
  for (e in x) {
    sum <- sum + (e*e);
  }
  sqrt(sum);
}

sampleHyperSphere <- function(dim) {
  x <- rnorm(dim, mean=0, sd=1);
  radius <- distance2(x);
  x/radius;
}

sampleHyperBall <- function(dim) {
  x <-  sampleHyperSphere(dim);
  x*runif(1);
}

sampleSimplexDirection <- function(dim) {
  x <- sampleHyperBall(dim-1);
  append(x, 0-sum(x));
}
