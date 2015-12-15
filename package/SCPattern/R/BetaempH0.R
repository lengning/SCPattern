#' @title Compare empirical non-directional KS statistics to KS statistics of permuted data
#' @usage BetaempH0(Data, DataShuffle,StateNames=c("EE","DE"),PiIn=NULL,Num=-100, iter=5, UpdatePi=TRUE)
#' @param Data directional KS statistics of empirical data
#' @param DataShuffle directional KS statistics of permuted data
#' @param PiIn initial value for Pi
#' @param Num The underflowed values will be imputed as exp(Num)
#' @param iter number of iterations
#' @param UpdatePi whether update Pi
#' @param StateNames names of the states
#' @examples tmp1 <- runif(10)
#' mat1 <- cbind(tmp1, 1-tmp1)
#' tmp2 <- runif(10)
#' mat2 <- cbind(tmp2, 1-tmp2)
#' res <- BetaempH0(mat1, mat2)
#' @author Ning Leng
#' @return Output: posterior probabilities and estimated hyper-parameters

BetaempH0 <- function(Data, DataShuffle,StateNames=c("EE","DE"),
 PiIn=NULL,Num=-100, iter=5, UpdatePi=TRUE){

res <- beta.mom(DataShuffle)
al <- res[1]
be <- res[2]
Mat <- rbind(c(al,be),c(be,al))


if(is.null(PiIn))PiIn <- rep(1/2,2)
P <- sapply(1:2,function(i){
	tmp <- dbeta(Data,Mat[i,1],Mat[i,2])
	tmp[which(is.na(tmp))] <- exp(Num)
tmp
})

Pi.trace <- PiIn
Pi <- PiIn

for(i in 1:iter){
PPi <- t(t(P)*Pi)
isNA <- which(is.na(PPi[,1]))
PP <- PPi/rowSums(PPi)
isNA <- which(is.na(PP[,1]))
if(length(isNA) > 0)PiNew <- colMeans(PP[-isNA,])
if(length(isNA) == 0)PiNew <- colMeans(PP)
if(UpdatePi==TRUE)Pi <- PiNew
Pi.trace <- rbind(Pi.trace,Pi)
}
rownames(P)=rownames(PP)=names(Data)
colnames(P)=colnames(PP)=StateNames
mm <- apply(PP,1,function(i)which.max(i)[1])
MAP <- colnames(PP)[mm]
maxPP <- sapply(1:length(mm),function(i)PP[i,mm[i]])
names(MAP)=names(maxPP)=names(Data)

Out <- list(PP=PP,P=P,Pi.trace=Pi.trace,
						al=al, be=be, par=Mat, MAP=MAP, maxPP=maxPP)
}

