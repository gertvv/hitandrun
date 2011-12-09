\name{har}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
"Hit and Run" sampler
}
\description{
The "Hit and Run" method generates a Markov Chain whose stable state converges on the uniform distribution over a convex polytope defined by a set of linear constraints.
}
\usage{
har(x0, constr, N, thin=1, homogeneous=FALSE, transform=NULL)
}
\arguments{
	\item{x0}{Starting point (must be in the polytope)}
	\item{constr}{Constraint definition (see details)}
	\item{N}{Number of iterations to run}
	\item{thin}{Thinning factor (keep every 'thin'-th sample)}
	\item{homogeneous}{Whether x0, constr and transform are given in homogeneous coordinate representation (see details)}
	\item{transform}{Transformation matrix to apply to the generated samples (optional)}
}
\details{
The constraints, starting point and transformation matrix can be given in homogeneous coordinate representation (an extra component is added to each vector, equal to 1.0).
This enables affine transformations (such as translation) to be applied to the coordinate vectors by the constraint and transformation matrices.
Be aware that while non-affine (perspective) transformations are also possible, they will not in general preserve uniformity of the generated samples.

Constraints are given as a list(constr=A, rhs=b, dir=d), where d should contain only "<=" (other constraints may be supported in the future).
The constraints define the polytope as usual for linear programming: Ax <= b.
In particular, it must be true that Ax0 <= b.
}
\value{
  A matrix containing the generated samples as rows.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Gert van Valkenhoef
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{generateSeedPoint}, \link{simplex.createTransform}, \link{simplex.createConstraints}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
