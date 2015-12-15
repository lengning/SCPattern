#' @title Compare empirical non-directional KS statistics to KS statistics of permuted data
#' @usage BetaempH0(Data, DataShuffle,PiIn=NULL,Num=-100, iter=5, UpdatePi=TRUE)
#' @param Data directional KS statistics of empirical data
#' @param DataShuffle directional KS statistics of permuted data
#' @param PiIn initial value for Pi
#' @param Num The underflowed values will be imputed as exp(Num)
#' @param iter number of iterations
#' @param UpdatePi whether update Pi
#' @param StateNames names of the states
#' @param one of DirReg, VGAM and optim. It defines the method to use to estimate hyper-parameters
#' @examples tmp1 <- runif(10)
#' mat1 <- cbind(tmp1, 1-tmp1)
#' tmp2 <- runif(10)
#' mat2 <- cbind(tmp2, 1-tmp2)
#' res <- BetaempH0(mat1, mat2)
#' @author Ning Leng
#' @return Output: posterior probabilities and estimated hyper-parameters


DirichletempH0 <- function(Data, DataShuffle,StateNames=c("EE","Down","Up"),
al=1, al2=1,be=10, PiIn=NULL, NumPat=3,Num=-100, iter=5, UpdatePi=T,
method="DirReg"){

DataShuffle0 <- DataShuffle
DataShuffle[which(Data==0)] <- exp(Num)

if(method=="optim"){
res <- optim(c(al,be),DirichletOptim, InputPool=list(DataShuffle),lower=c(0,0,0))
co <- c(res$par[1], res$par)
}
if(method=="VGAM"){
library(VGAM)
ydata <- data.frame(DataShuffle)
colnames(ydata) <- paste("y", 1:3, sep = "")
fit <- vglm(cbind(y1, y2, y3)  ~ 1, dirichlet, ydata, trace = TRUE, crit = "coef")
co <- exp(coef(fit))
}
if(method=="DirReg"){
#library(DirichletReg)
ydata <- data.frame(DataShuffle)
colnames(ydata) <- paste("y", 1:3, sep = "")
ALake <- data.frame(ydata)
ALake$Y <- DR_data(ydata[,1:3])
res1 <- DirichReg(Y ~ 1, ALake)
co <- exp(res1$coefficient)
}
al <- mean(co[1:2])
al2 <- mean(co[1:2])
be <- co[3]
message("alpha ", round(be,2), " beta ", round(al,2))
#library("gtools")
if(NumPat==4)
Mat <- rbind(c(al,al2,be), c(al,be,al2),c(be,al2,al),c(be,be,al))
if(NumPat==3)
Mat <- rbind(c(al,al2,be), c(al,be,al2),c(be,al2,al))
if(NumPat==2)
Mat <- rbind(c(al,be,al2),c(be,al2,al))

Data0 <- Data
Data[which(Data==0)] <- exp(Num) # 0 gives Inf
if(is.null(PiIn))PiIn <- rep(1/NumPat,NumPat)
P <- sapply(1:NumPat,function(i){
	tmp <- ddirichlet(Data,Mat[i,])
	tmp[which(is.na(tmp))] <- exp(Num)
tmp
})

rM <-rowMeans(P)

if("EE"%in%StateNames){
	P[which(rM==exp(Num)),] <- 0
	P[which(rM==exp(Num)),which(StateNames=="EE")] <- 1
}

Pi.trace <- PiIn
Pi <- PiIn

for(i in 1:iter){
PPi <- t(t(P)*Pi)
isNA <- which(is.na(PPi[,1]))
PP <- PPi/rowSums(PPi)
isNA <- which(is.na(PP[,1]))
if(length(isNA)>0)PiNew <- colMeans(PP[-isNA,])
if(length(isNA)==0)PiNew <- colMeans(PP)
if(UpdatePi==T)Pi <- PiNew
Pi.trace <- rbind(Pi.trace,Pi)
}
rownames(P)=rownames(PP)=rownames(Data)
colnames(P)=colnames(PP)=StateNames
mm <- apply(PP,1,function(i)which.max(i)[1])
MAP <- colnames(PP)[mm]
maxPP <- sapply(1:length(mm),function(i)PP[i,mm[i]])
names(MAP)=names(maxPP)=rownames(Data)
Out <- list(PP=PP,P=P,Pi.trace=Pi.trace,
alpha=be, beta=al, ab2=al2, par=Mat, MAP=MAP, maxPP=maxPP)
}

