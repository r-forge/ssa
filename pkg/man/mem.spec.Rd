\name{mem.spec}
\alias{mem.spec}
\title{ Maximum-entropy spectral estimate }
\description{
  Computes the maximum-entropy spectral estimate
}
\usage{ mem.spec(y,ord=5,...)
}
\arguments{
  \item{y}{ time series.}
  \item{ord}{ maximum order of the AR series to be estimated for y. }
  \item{...}{ Additional arguments passed to spec.ar}

}
\details{
  Uses spec.ar to compute the maximum-entropy estimate using Burg's method. See spec.ar
}
\value{
  
  \item{sp }{Description of 'comp1'}
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

  function(y, ord,... ) {
	
	y <- y - mean(y)
	T <- length(y)
	if ( ord >= T ) stop("The length of the window is longer than the data")

	sp<-spec.ar(y, method="burg",aic=FALSE, order=ord,...)
	return(sp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }