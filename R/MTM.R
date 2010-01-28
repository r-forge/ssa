
##	DPSS multitaper
##

mtm<-function(X,bdw=0.05,order=16,normalize=TRUE,taper=c("dpss")) {
  lx <- length(X)
  
  if (missing(taper))
    h<-taper(type="dpss",bandwidth=0.05,n.taper=8,n.sample=lx,normalize=TRUE)
  else
    h<-taper(type=taper,bandwidth=bdw,n.taper=order,n.sample=lx,normalize=TRUE)

#

Y<-h*matrix(rep(X,order),dim(h),byrow=T)
I<-t(abs(mvfft(t(Y)))^2)/(2*pi)
f<-apply(I,2,cumsum)
f<-f/matrix(rep(1:nrow(f),ncol(f)),dim(f))
  return(list(f=f,h=h,I=I))
}



### adaptive multitaper method
adp.mtm<-function(X,bdw=0.05,order=16,taper=c("dpss"),normalize=TRUE) {
  
  lx <- length(X)
  
  if (missing(taper))
    h<-taper(type="dpss",bandwidth=0.05,n.taper=8,n.sample=lx,normalize=TRUE)
  else
    h<-taper(type="dpss",bandwidth=bdw,n.taper=order,n.sample=lx,normalize=TRUE)

    
  Y<-h*matrix(rep(X,order),dim(h),byrow=T)
  I<-t(abs(mvfft(t(Y)))^2)/(2*pi)
  
  f.i<-apply(I[1:min(2,order),],2,mean)
  i<-outer(1:T,1:T,"-")
  
  A<-ifelse(i==0,2*L,sin(2*pi*L*i)/(pi*i))
  B<-matrix(rep(diag(h%*%A%*%t(h)),T),order,T)
  
  frel<-f.i/mean(f.i)
  frel<-matrix(rep(frel,order),order,T,byrow=T)
  B<-B/(B+(1-B)/frel)^2
  
  f<-apply(B*I,2,sum)/apply(B,2,sum)
  
  return(list(f=f,B=B))
}



##	Plotting calls...
##
plot.mtm <- function( ) {


}
