\name{simplex.sample}
\alias{simplex.sample}
\title{
Sample uniformly from a simplex
}
\description{
Generates uniform random variates over the (n-1)-simplex in n-dimensional space.
}
\usage{
simplex.sample(n, N)
}
\arguments{
	\item{n}{Dimension of the space}
	\item{N}{Number of samples to generate}
}
\details{
The samples will be uniform over the (n-1)-simplex.
}
\value{
  A matrix containing the generated samples as rows.
}
\references{
}
\author{
Gert van Valkenhoef
}
\note{
%%  ~~further notes~~
}
\seealso{
}
\examples{
n <- 3
N <- 10000
samples <- simplex.sample(n, N)

# Check dimension
stopifnot(dim(samples) == c(N, n))

# Check that w_i >= 0
stopifnot(samples >= 0)

# Check that sum_i w_i = 1
E <- 1E-12
stopifnot(apply(samples, 1, sum) > 1 - E)
stopifnot(apply(samples, 1, sum) < 1 + E)
}
\keyword{rejection}
\keyword{uniform}
\keyword{sampling}
