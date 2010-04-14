\name{btpsd}
\alias{btpsd}

\title{ Blackman-Tukey spectral density estimate }
\description{
	Computes the Blackman-Tukey spectral density estimate}
\usage{
}
\arguments{
  \item{y}{ time series. }
  \item{type}{ Character. Default is "Tukey". Other window choices are "Triangular","Hanning",and "Hamming".}
  \item{win}{ window length.}
  \item{taper}{ amount to taper the time series. Default is 0.5.}

}
\details{
}
\value{
  \item{est }{ Blackman-Tukey estimates of the power spectral density.}
}
\references{ D.B. Percival and A. Walden (1993), Spectral Analysis for Physical Applications: Multitaper and Conventional Univariate Techniques, Cambridge University Press, Cambridge, UK.}
\author{ Patrick Crutcher}
\note{ ~~further notes~~ 

}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(y, type="Tukey", win=N , taper=0.5) {
  
  y <- y - mean(y)
  T <- length(y)
  if ( win >= T ) stop("The length of the window is longer than the data length.")

  if (is.null(win)) N <- ceiling(2*sqrt(T)) else N <- win
  
  if (type =="Tukey") { w <- tukey(win,taper) }
  if (type =="Hanning") { w <- hanning.window(win) }
  if (type =="Hamming") { w <- hamming.window(win) }
  if (type =="Triangular") { w <- tri.window(win) }


	r <- acf(y, lag.max = N-1,plot=FALSE)$acf;

	rw <- r*w[-N]
	
	est <- Mod(fft(rw))^2/(2*pi*length(y))
	return(est)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }