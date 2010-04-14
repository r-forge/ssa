\name{mtm}
\alias{mtm}

\title{ Multi-taper method }
\description{
	Computes the multi-taper estimate of the power spectral density}
\usage{
mtm(y, bdw = 0.05, order = 10, normalize = TRUE, tape =
                 c("dpss"))}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ a vector containing a time series of data }
  \item{bdw}{ bandwidth for the Diploid Prolate Spherical Sequences.}
  \item{order}{ an integer defining the number of tapers to use with the Prolate sequences.}
  \item{normalize}{ Logical. If TRUE, taper will be normalized. }
  \item{tape}{ a character string denoting the type of taper to create. See \pkg{sapa}. Default: "dpss"}

}
\details{
  Utilizes the \pkg{sapa} library for tapers.
  }
\value{

  \item{f}{vector of spectral density function estimates}
  \item{h }{vector of tapers}
  \item{I}{ periodogram estimates for each segment, from 1 to order.}
  
  }
\references{ D.B. Percival and A. Walden (1993), Spectral Analysis for Physical Applications: Multitaper and Conventional Univariate Techniques, Cambridge University Press, Cambridge, UK.}
\author{ Patrick Crutcher }
\note{ ~~further notes~~ 
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
  #		y <- data(tsoi)
  #      mtm(y,bdw=0.05,order=10,normalize=TRUE,tape=c("dpss"))
  #      mtm(y,bdw=0.025,order=5,normalize=FALSE,tape=c("sine"))
  
## The function is currently defined as
function(y,bdw=0.05,order=10,normalize=TRUE,tape=c("dpss")) {

  ly <- length(y)
  
  if (missing(tape))
    h<-taper(type="dpss",bandwidth=0.05,n.taper=8,n.sample=ly,normalize=TRUE)
  else
    h<-taper(type=tape,bandwidth=bdw,n.taper=order,n.sample=ly,normalize=TRUE)

#

Y<-as.matrix(h)*matrix(rep(y,order),dim(h),byrow=TRUE)
I<-t(abs(mvfft(t(Y)))^2)/(2*pi)
f<-apply(I,2,cumsum)
f<-f/matrix(rep(1:nrow(f),ncol(f)),dim(f))
  return(list(f=f,h=h,I=I))
}
\keyword{ ts }