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
E <- 1E-15 # Error limit

sample <- hypersphere.sample(n, N)

# Check summing to unity
all(abs(apply(hypersphere.sample(3, 2), 1, function(x) {sum(x^2) })) - 1
< E)
}
\keyword{hypersphere}
\keyword{uniform}
\keyword{sampling}
