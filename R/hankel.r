##Massimo Bilancia
##22092008

function(X){

if(!is.matrix(X))
stop("I need a matrix!")
L <- dim(X)[1]
K <- dim(X)[2]
hankel.matrix <- matrix(0, nrow=L, ncol=K)
shift.block <- matrix(0, nrow=L, ncol=K+(L-1))

for(i in 1:L)
{ shift.block[i,i:((K-1)+i)] <- X[i,] }	
nonzero.counts <- function(x){ l <- length(x[x!=0]) ; return(l)}
temp.sum <- apply(shift.block, 2, sum)
temp.length <- apply(shift.block, 2, nonzero.counts)
average <- temp.sum/temp.length

for(i in 1:(K+(L-1)))
{ shift.block[,i] <- average[i] }

for(i in 1:L)
{ hankel.matrix[i,]<- shift.block[i,i:(K+(i-1))] }

return(hankel.matrix)
}
