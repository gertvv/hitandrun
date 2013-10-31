\name{hypersphere.sample}
\alias{hypersphere.sample}
\title{
Sample uniformly from an n-hypersphere
}
\description{
Generates uniform random variates over an n-hypersphere
}
\usage{
hypersphere.sample(n, N)
}
\arguments{
  \item{n}{Dimension of the hypersphere}
  \item{N}{Number of samples}
}
\value{
	A single n-dimensional sample from the hypersphere.
}
\author{
Tommi Tervonen <tommi@smaa.fi>
}
\examples{
n <- 3 # Dimension
N <- 5 # Nr samples

sample <- hypersphere.sample(n, N)

# Check summing to unity
vec.norm <- function(x) { sum(x^2) }
stopifnot(all.equal(apply(sample, 1, vec.norm), rep(1, N)))
}
\keyword{hypersphere}
\keyword{uniform sampling}
