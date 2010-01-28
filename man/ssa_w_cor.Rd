\name{ssa_w_cor}
\alias{ssa_w_cor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ssa_w_cor(z, l, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{ ~~Describe \code{z} here~~ }
  \item{l}{ ~~Describe \code{l} here~~ }
  \item{k}{ ~~Describe \code{k} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
