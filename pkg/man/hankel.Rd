\name{hankel}
\alias{hankel}

\title{ Apply Hankel operator to a time series }
\description{
	Computes trajectory matrix as a Hankel matrix.}
\usage{
hankel(x, L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ time series }
  \item{L}{ window length}
}
\details{
  Calls C code to compute the Hankel matrix
}
\value{

  \item{value }{a matrix in Hankel form}
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
function(x,L){
cc<-rep(0,L*L); T<-length(x)
cl<-.C("hankelC",
        as.double(x),
        as.double(cc),
        as.integer(L),
        as.integer(T))
return(matrix(cl[[2]],L,L))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }