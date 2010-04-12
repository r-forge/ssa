
mtm<-function(y,bdw=0.05,order=10,normalize=TRUE,tape=c("dpss")) {

  # Description:
  #   Performs multi-taper method on a time series
  #
  # Usage:
  #
  #   mtm(y,bdw=0.05,order=10,normalize=TRUE,taper=c("dpss"))
  #
  # Arguments:
  #
  #   y:      a vector containing a time series of data
  #
  #   bdw:    bandwidth for the Diploid Prolate Spherical Sequences. 
  #
  #   order:	an integer defining the number of tapers to use with the Prolate sequences. 
  #
  #   normalize:    Logical. If TRUE, taper will be normalized. 
  #
  #	  tape: 	a character string denoting the type of taper to create. See library(sapa). Default: "dpss"
  #
  # Value:
  # 
  #   f:	vector of spectral density function estimates
  #
  #	  h: 	vector of tapers
  #
  #	  I: 	periodogram estimates for each segment, k from 1 to order.
  #
  # Author:
  #
  #   Patrick Crutcher
  # 
  # Examples:
  #
  #      mtm(y,bdw=0.05,order=10,normalize=TRUE,tape=c("dpss"))
  #      mtm(y,bdw=0.025,order=5,normalize=FALSE,tape=c("sine"))
  #
  
  lx <- length(y)
  
  if (missing(tape))
    h<-taper(type="dpss",bandwidth=0.05,n.taper=8,n.sample=lx,normalize=TRUE)
  else
    h<-taper(type=tape,bandwidth=bdw,n.taper=order,n.sample=lx,normalize=TRUE)

#

Y<-as.matrix(h)*matrix(rep(y,order),dim(h),byrow=T)
I<-t(abs(mvfft(t(Y)))^2)/(2*pi)
f<-apply(I,2,cumsum)
f<-f/matrix(rep(1:nrow(f),ncol(f)),dim(f))
  return(list(f=f,h=h,I=I))
}



### adaptive multitaper method
adp.mtm<-function(y,bdw=0.05,order=10,normalize=TRUE,tape=c("dpss")) {
  
  # Description:
  #   Performs Adaptive multi-taper method on a time series
  #
  # Usage:
  #
  #   adp.mtm(y,bdw=0.05,order=10,normalize=TRUE,taper=c("dpss"))
  #
  # Arguments:
  #
  #   y:      a vector containing a time series of data
  #
  #   bdw:    bandwidth for the Diploid Prolate Spherical Sequences. 
  #
  #   order:	an integer defining the number of tapers to use with the Prolate sequences. 
  #
  #   normalize:    Logical. If TRUE, taper will be normalized. 
  #
  #	  tape: 	a character string denoting the type of taper to create. See library(sapa). Default: "dpss"
  #
  # Value:
  # 
  #   f:	vector of spectral density function estimates
  #
  #	  B:  	weighting function values
  #
  #   I: 	periodogram estimates for each segment, k from 1 to order.
  #
  # Author:
  #
  #   Patrick Crutcher
  # 
  # Examples:
  #
  #      mtm(y,bdw=0.05,order=10,normalize=TRUE,tape=c("dpss"))
  #      mtm(y,bdw=0.025,order=5,normalize=FALSE,tape=c("sine"))
  #
  
  ly <- length(y)
  
  if (missing(tape))
    h<-taper(type="dpss",bandwidth=0.05,n.taper=8,n.sample=ly,normalize=TRUE)
  else
    h<-taper(type=tape,bandwidth=bdw,n.taper=order,n.sample=ly,normalize=TRUE)

  hh<-as.matrix(h)
    
  Y<- hh*matrix(rep(y,order),dim(hh),byrow=TRUE)
  I<-t(abs(mvfft(t(Y)))^2)/(2*pi)
  
  f.i<-apply(I[1:min(2,order),],2,mean)
  i<-outer(1:ly,1:ly,"-")
  
  A<-ifelse(i==0,2*bdw,sin(2*pi*bdw*i)/(pi*i))
  B<-matrix(rep(diag(hh%*%A%*%t(hh)),ly),order,ly)
  
  frel<-f.i/mean(f.i)
  frel<-matrix(rep(frel,order),order,ly,byrow=TRUE)
  B<-B/(B+(1-B)/frel)^2
  
  f<-apply(B*I,2,sum)/apply(B,2,sum)
  
  return(list(f=f,B=B,I=I))
}
