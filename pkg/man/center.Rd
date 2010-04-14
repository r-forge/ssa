\name{center}
\alias{center}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Centering of a series }
\description{ Centers a time series
}
\usage{
center(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ time series or vector. }
}
\details{ Returns the series minus the mean.
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{value }{returns the centered series.}
}
\references{ ~put references to the literature/web site here ~ }
\author{ Patrick Crutcher }
\note{ ~~further notes~~ 

}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x) x-mean(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }