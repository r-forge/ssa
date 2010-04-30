detrend <- function(y,demean=TRUE,det=FALSE,norm=FALSE, smooth=FALSE,...){

  # Description:
  #   Preprocessors for a time series
  #
  # Usage:
  #
  #   detrend(y,demean=TRUE,det=FALSE,norm=FALSE,smooth=FALSE,...)
  #
  # Arguments:
  #
  #   y:      a vector containing a time series of data
  #
  #   demean: Logical.If TRUE, time series y will be centered. 
  #
  #   det: 	  Logical.If TRUE, then time trend is removed by fitting a local weighted
  #           regression (loess) to the time series y
  #
  #   smooth: Logical.If TRUE, then time trend is removed by a local fitting regression(loess) 
  #			  to the time series y
  #
  #   ...:    Additional arguments passed to loess
  #
  # Value:
  # 
  #   Vector with the detrended data.
  #
  # Author:
  #
  #   Patrick Crutcher
  # 
  # Examples:
  #
  #      ts <- rnorm(100,25,2)
  #      detrend(ts)
  #      detrend(ts,smooth=TRUE)
  #      detrend(ts,smooth=TRUE,f=1/3)
	


	if(center==TRUE){
		
		dt<- y - mean(y)
		return(dt)
	}	
	
	if(norm==TRUE){
		
		dt<- (y - mean(y))/var(y)
		return(dt)
	}

	if(det==TRUE){
		
		x<- seq(1:length(y))
		dt<- lm(y ~ x)$res
		return(dt)
		}
	if(smooth==TRUE){
	
		dt<- y-loess(y,...)$y
		return(dt)
	}
}