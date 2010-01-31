
btpsd <- function(y, type="Tukey", win=N , taper=0.5) {
	## INPUT ##
	## y = time series of interest
	## win = window length = number of autocorrelations in estimation, DEFAULT = 2 * sqrt(lenth(y))
    ## default window =  Blackman-Tukey
  
  y <- y - mean(y)
  T <- length(y)
  if ( win >= T ) stop("The length of the window is longer than the data length.")

  if (is.null(win)) N <- ceiling(2*sqrt(T)) else N <- win
  
  if (type =="Tukey") { w <- tukey(win,taper) }
  if (type =="Hanning") { w <- hanning.window(win) }
  if (type =="Hamming") { w <- hamming.window(win) }
  if (type =="Triangular") { w <- tri.window(win) }


	r <- acf(y, lag.max = N-1,plot=FALSE)$acf;

	rw <- r*w[-N]
	
	est <- Mod(fft(rw))^2/(2*pi*length(y))
	return(est)
}

mem.spec <- function(y, ord,... ) {
	## INPUT ##
	## y = time series of interest
	## order of AR series
	## ... any params to pass to spec.ar
  
	y <- y - mean(y)
	T <- length(y)
	if ( ord >= T ) stop("The length of the window is longer than the data")

	sp<-spec.ar(y, method="burg",aic=FALSE, order=ord,...)
	return(sp)
}

tukey<-function(n,a) 
{   
	t2 <- seq(0,1,length.out=n);
    per <- a/2; 
    tl <- floor(per*(n-1))+1;
    th <- n-tl+1;
    # Window is defined in three sections: taper, constant, taper
    w <- c( 0.5+0.5*cos( (pi/per)*(t2[1:tl] - per) ) ,  rep(1,th-tl-1), 0.5+0.5*cos((pi/per)*(t2[th:n] - 1 + per))) 
	return(w)
}

hanning.window<-function(n) 
{
    if (n == 1 ) 
        c <- 1
    else {
        n <- n - 1
        c <- 0.5 - 0.5 * cos(2 * pi * (0:n)/n)
    }
    return(c)
}

tri.window<-function(n)
{
	c <- c( 2/n*( n/2-abs( seq(0,n-1)-(n-1)/2 ) ) )

	return(c)
}



hamming.window<-function(n) 
{
    if (n == 1) 
        c <- 1
    else {
        n <- n - 1
        c <- 0.54 - 0.46 * cos(2 * pi * (0:n)/n)
    }
    return(c)
}
