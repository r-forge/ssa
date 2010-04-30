#
#   ssa package
#   Copyright (C) 2008  Jan de Leeuw <deleeuw@stat.ucla.edu>
#   UCLA Department of Statistics, Box 951554, 
#   Los Angeles, CA 90095-1554
#   
#   This program is free software; you can redistribute it 
#   and/or modify it under the terms of the GNU General Public 
#   License as published by the Free Software Foundation; 
#   either version 2 of the License, or (at your option) 
#   any later version.
#
#   This program is distributed in the hope that it will be 
#   useful, but WITHOUT ANY WARRANTY; without even the implied 
#   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
#   PURPOSE.  See the GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public 
#   License along with this program; if not, write to the 
#   Free Software Foundation, Inc., 675 Mass Ave, Cambridge, 
#   MA 02139, USA.
#
###################################################################
#
# version 0.0.1, 2008-10-10, initial
# version 0.1.0, 2008-10-16, similarity and cluster parameters
# version 0.2.0, 2008-10-17, got rid of hankel
# version 0.3.0, 2008-10-17, added Toeplitz SSA
# version 0.3.1, 2008-10-19, C code for crossproduct
# version 0.3.2, 2008-10-19, preprocessing option
# version 0.3.3, 2008-10-19, C code for fromCross
# version 0.4.0, 2008-10-19, added Fourier SSA
#
# to do: 
#
# -- add clustering methods
# -- add similarity measures
# -- add preprocessors
# -- ndim for eigen
# -- better algorithm for Toeplitz eigenvectors
#

.First.lib <- function(lib, pkg) {
  library.dynam("ssa", pkg, lib)
}

ssa<-function(data,L=floor(length(data)/2),b=as.list(1:L),preproc=c(),method="hankel",ssa_similarity=ssa_w_cor,ssa_cluster=trans,par=.25,...) {
    if (!identical(class(data),"ts")) stop("Data should be of class ts")
    
    if ( !identical(preproc,"NULL")) {   
        data<-detrend(data,preproc,...);
        }
    
    T<-length(data); K<-T-L+1
    
    if (identical(method,"toeplitz")) {   
        cross<-toeplitz(rconv(data,1:L))
        s<-eigen(cross)$vectors
        }
    if (identical(method,"hankel")) {
        cross<-hankel(data,L)
        s<-eigen(cross)$vectors
        }
    if (identical(method,"fourier")) {
        s<-sincos(L)
        }
    a<-fromCross(data,s)
    
    r<-ssa_similarity(a,L,K)
    
    if (identical(ssa_cluster,"trans")) {
        b<-ssa_trans(a,r)
        }
    if (identical(ssa_cluster,"hierch")) {
        b<-ssa_hclust(r,par)
        }
   # b<-ssa_cluster(r,par)
    
    a<-sapply(b,function(kk) rowSums(as.matrix(a[,kk])))
    y<-ts(a,start=start(data),end=end(data),frequency=frequency(data))
    return(list(y=y,b=b,s=colSums(a^2),r=r))
}

ssa_w_cor<-function(z,l,k) {
    ls<-min(l,k); ks<-max(l,k)
    n<-nrow(z); w<-rep(ls,n)
    w[1:ls]<-1:ls; w[(ks+1):n]<-n-(ks:(n-1))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
    c<-crossprod(z,w*z); d<-diag(c)
    return(c/sqrt(outer(d,d)))
}

### needs work ###
ssa_hclust<-function(r,par,a,L,K) {

	if (identical(r,"NULL")) {
        r<-as.dist(1-abs(ssa_w_cor(a,L,K)))
        }
	s<-hclust( as.dist(r), method = par[1])
	return(s)

}

ssa_trans<-function(r,cut) {
    s<-ifelse(abs(r)>cut,1,0)
    v<-warshall(s)
    h<-unique(v)
    return(apply(h,1,function(v) which(v==1)))
}

warshall<-function(a) {
n<-nrow(a)
for (j in 1:n) {
    for(i in 1:n) {
        if (a[i,j]==1) a[i,]<-pmax(a[i,],a[j,])
        }
    }
return(a)
}

sincos<-function(L) {
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

rconv<-function(x,lag) {
.C("cconv",
    as.integer(lag),
    as.integer(length(lag)),
    as.double(x),
    as.integer(length(x)),
    as.double(vector("double",length(lag))))[[5]]
}

hankel<-function(x,L){
cc<-rep(0,L*L); T<-length(x)
cl<-.C("hankelC",
        as.double(x),
        as.double(cc),
        as.integer(L),
        as.integer(T))
return(matrix(cl[[2]],L,L))
}

fromCross<-function(x,s) {
    L<-nrow(s); T<-length(x)
    aa<-rep(0,T*L); s<-as.vector(s)
    al<-.C("fromCrossC",
            as.double(x),
            as.double(s),
            as.double(aa),
            as.integer(L),
            as.integer(T))
    return(matrix(al[[3]],T,L))
}

is.even<-function(x) return((as.integer(x) %% 2) == 0)

center<-function(x) x-mean(x)

ident<-function(x) x

