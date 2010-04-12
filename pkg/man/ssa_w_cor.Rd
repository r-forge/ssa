\name{ssa_w_cor}
 \alias{ssa_w_cor}

\title{ w-correlation }
\description{
  Computes w - correlations (Golyandina, et.al. 2001) from a matrix
  containing reconstructed components in its columns.
}
\usage{
ssa_w_cor(z, l, k)
}
\arguments{
  \item{z}{ A matrix containing reconstructed components in its columns. }
  \item{l}{ An integer for the window width }
  \item{k}{ An integer defined by k=t-l+1 }
}
\details{
  w-correlations can be used to assess how well reconstructed components can be separated from each other. 
  See reference for more details.
}
\value{
    A square matrix containing the w-correlations between components.
}
\references{ Golyandina, N.; Nekrutkin, V. & Zhiglkilavskifi, A. Analysis of Time Series Structure: SSA and Related Techniques. CRC Press, 2001
 }
\author{ Patrick Crutcher }
\note{ ~~further notes~~ 
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

#x <- sin(seq(0,10*pi,len=200))
#x <- x + rnorm(x)/2
#x.wc <- w.cor(x.rc,40)
#image(x.wc,col=gray(100:0/100))

## The function is currently defined as
function(z,l,k) {
    ls<-min(l,k); ks<-max(l,k)
    n<-nrow(z); w<-rep(ls,n)
    w[1:ls]<-1:ls; w[(ks+1):n]<-n-(ks:(n-1))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
    c<-crossprod(z,w*z); d<-diag(c)
    return(c/sqrt(outer(d,d)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }