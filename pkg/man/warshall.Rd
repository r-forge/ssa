\name{warshall}
\alias{warshall}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Warshall's algorithm }
\description{
  Warshall's algorithm to compute the transitive closure for an adjacency matrix.
}
\usage{
warshall(a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{ adjacency matrix }
}
\details{

}
\value{
 
  \item{a}{ transitive closure of the relation corresponding to a }
}
\references{ S. Warshal, A Theorem on Boolean Matrices. Journal of the Association of Computer Machinery, 9:11-12,1962}
\author{ Patrick Crutcher }
\note{ ~~further notes~~ 

}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(a) {
n<-nrow(a)
for (j in 1:n) {
    for(i in 1:n) {
        if (a[i,j]==1) a[i,]<-pmax(a[i,],a[j,])
        }
    }
return(a)
}
\keyword{warshall}