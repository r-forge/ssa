\name{sincos}
\alias{sincos}
\title{ Fourier SSA }
\description{
  Fourier SSA
}
\usage{
sincos(L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{ window length }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{s}{matrix of Fourier coefficients}

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
function(L) {
    x<-2*pi*(1:L)/L
    M<-floor(L/2)
    if (is.even(L)) N<-M-1 else N<-M
    s<-matrix(0,L,N+M+1)
    for (k in 1:N)
        s[,k]<-sin(x*k)
    for (k in 0:M)	        
        s[,N+k+1]<-cos(x*k)
return(apply(s,2,function(z) z/sqrt(sum(z^2))))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}