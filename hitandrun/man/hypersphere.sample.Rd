\name{hypersphere.sample}
\alias{hypersphere.sample}
\title{
Sample uniformly from an n-hypersphere
}
\description{
Generates uniform random variates over an n-hypersphere
}
\usage{
hypersphere.sample(n)
}
\arguments{
	\item{n}{Dimension of the hypersphere}
}
\value{
	A single n-dimensional sample from the hypersphere.
}
\author{
Tommi Tervonen <tommi@smaa.fi>
}
\examples{
n <- 3
E <- 1E-15 # Error limit

sample <- hypersphere.sample(3)

# Check summing to unity
stopifnot(abs(sum(sample^2) - 1) < E)

}
\keyword{hypersphere}
\keyword{uniform}
\keyword{sampling}
