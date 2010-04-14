\name{adp.mtm}
\alias{adp.mtm}

\title{ Performs Adaptive multi-taper method on a time series }
\description{
	Computes the adaptive multi-taper estimate of the power spectral density}
\usage{
adp.mtm(y,bdw=0.05,order=10,normalize=TRUE,tape=c("dpss"))
}
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
  \item{B}{weighting function values for the adaptive estimate}
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

## The function is currently defined as
library(sapa)
function(y,bdw=0.05,order=10,normalize=TRUE,tape=c("dpss")){
  
  ly <- length(y)
  
  if (missing(tape))
    h<-taper(type="dpss",bandwidth=0.05,n.taper=8,n.sample=ly,normalize=TRUE)
  else
    h<-taper(type=tape,bandwidth=bdw,n.taper=order,n.sample=ly,normalize=TRUE)

  hh<-as.matrix(h)
    
  Y<-hh*matrix(rep(y,order),dim(hh),byrow=TRUE)
  I<-t(abs(mvfft(t(Y)))^2)/(2*pi)
  
  f.i<-apply(I[1:min(2,order),],2,mean)
  i<-outer(1:ly,1:ly,"-")
  
  A<-ifelse(i==0,2*bdw,sin(2*pi*bdw*i)/(pi*i))
  B<-matrix(rep(diag(hh%*%A%*%t(hh)),ly),order,ly)
  frel <- f.i/mean(f.i)
  frel <- matrix(rep(frel,order),order,ly,byrow=TRUE)
  B<-B/(B+(1-B)/frel)^2
  
  f<-apply(B*I,2,sum)/apply(B,2,sum)
  
  return(list(f=f,B=B,I=I))
}

\keyword{ ts }