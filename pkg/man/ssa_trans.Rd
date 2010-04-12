\name{ssa_trans}
\alias{ssa_trans}
\title{  Equivalence classes for the transitive closure of a matrix}
\description{
   Equivalence classes under the transitive closure}

\usage{
ssa_trans(r, cut)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{ an adjancency matrix}
  \item{cut}{ threshold value }
}
\details{

}
\value{

  \item{ value }{equivalence classes used to define groups}
}
\references{  }
\author{ Patrick Crutcher }
\note{ ~~further notes~~ 

}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(r,cut) {
    s<-ifelse(abs(r)>cut,1,0)
    v<-warshall(s)
    h<-unique(v)
    return(apply(h,1,function(v) which(v==1)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }