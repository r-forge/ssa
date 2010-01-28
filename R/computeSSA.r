##Massimo Bilancia
##220908
function(x, L=60){

N <- length(x)
K <- N - L + 1

X <- matrix(0, nrow=L, ncol=K)
for(i in 1:L){
	for(j in 1:K){
		X[i,j]=x[i+j-1]}}

X.svd <- svd(X/sqrt(K))

ifelse(L<=K, sup.index<-L, sup.index<-K)
X.elem.mat <- array(data=0, dim=c(L,K,sup.index)) 
X.elem.comp.array <- array(data=0, dim=c(L,K,sup.index)) 
for(i in 1:sup.index){
X.elem.mat[,,i] <- (sqrt(K)*X.svd$d[i]) * (X.svd$u[,i] %*% t(X.svd$v[,i]))
}


X.elem.comp <- matrix(data=0, nrow=N, ncol=sup.index)
for(i in 1:sup.index){
temp <- hankel(X.elem.mat[,,i])
X.elem.comp.array[,,i] <- temp
series <- rbind(as.matrix(temp[,1]), as.matrix(temp[L,2:K]))
X.elem.comp[,i] <- series
}

return(list(SVD=X.svd, ARRAY=X.elem.comp.array, COMP=X.elem.comp))
}
