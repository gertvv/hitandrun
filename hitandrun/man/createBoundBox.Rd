\name{createBoundBox}
\alias{createBoundBox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate bounding box based on given constraint matrix.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
createBoundBox(constr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{constr}{Constraint definition
%%     ~~Describe \code{constr} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (constr) 
{
    n <- dim(constr$constr)[2]
    extreme <- findExtremePoints(constr)
    lb <- apply(extreme, 1, min)
    ub <- apply(extreme, 1, max)
    boundFn <- function(x, d) {
        c(max(sapply(1:(n - 1), function(i) {
            min((lb[i] - x[i])/d[i], (ub[i] - x[i])/d[i])
        })), min(sapply(1:(n - 1), function(i) {
            max((lb[i] - x[i])/d[i], (ub[i] - x[i])/d[i])
        })))
    }
    start <- (1/(2 * (n - 1))) * apply(extreme, 1, sum)
    list(bound = boundFn, start = start, lb = lb, ub = ub)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
