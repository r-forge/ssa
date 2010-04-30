\name{detrend}
\alias{detrend}

\title{ Detrend a time series }
\description{
	Preprocessors for a time series.}
\usage{
detrend(y,demean=TRUE,det=FALSE,norm=FALSE,smooth=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ time series }
  \item{demean}{ Logical.If TRUE, time series y will be centered. }
  \item{norm}{ Logical.If TRUE, time series y will be centered and normalized }
  \item{det}{ Logical.If TRUE, then time trend is removed by fitting a linear model to the time series y }
  \item{smooth}{ Logical.If TRUE, then time trend is removed by a local fitting regression(loess) to the time series y}
  \item{...}{ Additional arguments passed to loess}

}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{

  \item{dt }{Vector with the detrended data.}
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }