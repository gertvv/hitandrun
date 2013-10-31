\name{simplex.sample}
\alias{simplex.sample}
\title{
Sample uniformly from a simplex
}
\description{
Generates uniform random variates over the (n-1)-simplex in n-dimensional space.
}
\usage{
simplex.sample(n, N, sort=FALSE)
}
\arguments{
	\item{n}{Dimension of the space}
	\item{N}{Number of samples to generate}
	\item{sort}{Whether to sort the components in descending order}
}
\details{
The samples will be uniform over the (n-1)-simplex.
}
\value{
	\item{samples}{A matrix containing the generated samples as rows.}
}
\author{
Gert van Valkenhoef
}
\examples{
n <- 3
N <- 10000
samples <- simplex.sample(n, N)$samples

# Check dimension
stopifnot(dim(samples) == c(N, n))

# Check that w_i >= 0
stopifnot(samples >= 0)

# Check that sum_i w_i = 1
E <- 1E-12
stopifnot(apply(samples, 1, sum) > 1 - E)
stopifnot(apply(samples, 1, sum) < 1 + E)


## Now with descending order
samples <- simplex.sample(n, N, sort=TRUE)$samples

# Check dimension
stopifnot(dim(samples) == c(N, n))

# Check that w_i >= 0
stopifnot(samples >= 0)

# Check that sum_i w_i = 1
E <- 1E-12
stopifnot(apply(samples, 1, sum) > 1 - E)
stopifnot(apply(samples, 1, sum) < 1 + E)

# Check w_i >= w_{i+1}
stopifnot(samples[,1] >= samples[,2])
stopifnot(samples[,2] >= samples[,3])
}
\keyword{simplex}
\keyword{uniform sampling}
