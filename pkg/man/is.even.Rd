\name{is.even}
\alias{is.even}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Divisibility by 2   }
\description{
  Determines if an integer is even or not
}
\usage{
is.even(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ integer }
}
\details{
}
\value{
  
  \item{comp1 }{logical. If x is divisible by two, it returns TRUE. }
  ...
}
\references{  }
\author{ Patrick Crutcher}
\note{ ~~further notes~~ 

}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x) return((as.integer(x) \%\% 2) == 0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }