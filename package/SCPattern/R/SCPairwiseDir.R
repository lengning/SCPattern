#' @title Perform directional tests across two or more conditions
#' @usage SCPairwiseDir <- function(Data, Conditions, sizeFactors,
#' Dropout.remove=FALSE, Dropout.upper=0, maxround=5, NumPat=3,
#' StateNames=c("EE","Down","Up"),UpdatePi=TRUE, Print=TRUE,SmallNum=1,  Seed=10, LOD=50)
#' @param Data A data matrix contains expression values for each transcript
#'    (gene or isoform level). In which rows should be transcripts
#'    and columns should be samples. The matrix is expected to be unnormalized (raw data).
#' @param Conditions A factor indicates the condition which each sample belongs to.
#' @param sizeFactors The normalization factors. It should be a vector with lane
#'    specific numbers (the length of the vector should be the same
#'    as the number of samples, with the same order as the columns of Data).
#' @param Droupout.remove,Dropout.upper For a given gene, whether cells with dropout will be 
#' ignored in the test. Defaulte is FALSE. If it is set to TRUE, values that are less or equal to
#' Dropout.upper will be ignored during the test.
#' @param maxround Number of iterations. The default value is 5.
#' @param NumPat number of patterns to consider. Could be 3 or 2. Default is 3.
#' When NumPat = 3, three states (no change, down and up) are considered for each 
#' transition. When NumPat = 2, two stages (down and up) are considered.
#' When NumPat = 4, four stages (no change, down, up and both directions) are considered.
#' @param StateNames names of the states
#' @param UpdatePi default is TRUE
#' @param Print whether print the messages
#' @param SmallNum When calculating log2 expression, a small number was added to each of
#' the values. Default is 1.
#' @param Seed seed
#' @param LOD Genes with max value < LOD will be eliminated in the analysis.
#' @param method one of DirReg, VGAM and optim. It defines the method to use to estimate hyper-parameters
#' @examples mat <- matrix(exp(rnorm(100,10,1)),ncol=10)
#' Sizes <- MedianNorm(mat)
#' tmp <- SCPairwiseDir(mat, rep(c("c1","c2"),each=5), Sizes, LOD=0)
#' @author Ning Leng
#' @return Output: For each comparison, output posterior probabilities of being each state.

SCPairwiseDir <- function(
		Data, Conditions, sizeFactors,
		Dropout.remove=FALSE, Dropout.upper=0,
		maxround=5,
		NumPat=3,
		StateNames=c("EE","Down","Up"),
		UpdatePi=TRUE, Print=TRUE, 
		SmallNum=1, Seed=10, LOD=50, method="DirReg"
	){

set.seed(Seed)
Dataraw <- Data
DataNorm <- GetNormalizedMat(Data, sizeFactors)

Gt10 <- which(apply(DataNorm,1,max)>LOD)

NonAllZeroNames <- Gt10 
if(length(NonAllZeroNames)<nrow(Data) & Print==TRUE) cat("Removing low expressed transcripts \n")
if(length(NonAllZeroNames)>0) Data <- Data[NonAllZeroNames,]
if(!is.factor(Conditions))Conditions <- as.factor(Conditions)

	
NumCond <- nlevels(Conditions)
CondLevels <- levels(Conditions)

	

# Divide by SampleSize factor
DataList.unlist=DataListIn.unlist=Data
if(length(sizeFactors)==ncol(Data))
DataList.unlist.dvd <- t(t( DataList.unlist)/sizeFactors)

if(length(sizeFactors)!=ncol(Data))
DataList.unlist.dvd <- DataList.unlist/sizeFactors

DataList.unlist.dvd.log=log2(DataList.unlist.dvd + SmallNum)

###################################
# Calculate KS stat for all comparisons
###################################
TestRes <- KSDir(DataList.unlist.dvd.log, Conditions,Dropout.remove=Dropout.remove, 
			 Dropout.upper=log2(Dropout.upper+SmallNum))
List <- TestRes$List


# shuffle
Which <- sapply(1:NumCond,function(i)which(Conditions==CondLevels[i]),simplify=FALSE)
Shuffle <- sapply(Which,function(i)sample(i,length(i)))
ShuffleV <- unlist(Shuffle)
DataShuffle <- DataList.unlist.dvd.log[,ShuffleV]

ListRandom <- KSDir(DataShuffle, Conditions,Dropout.remove=Dropout.remove, 
	 Dropout.upper=Dropout.upper)$List


#################################
# impute 0s

List3D=ListRandom3D=vector("list",length(List))
for(i in 1:length(List)){
for(j in 1:2){
tmp1 <- ListRandom[[i]][,j]
Min <- min(tmp1[which(tmp1>0)])
ListRandom[[i]][which(tmp1<=0),j] <- Min

tmp2 <- List[[i]][,j]
Min2 <- min(tmp2[which(tmp2>0)])
List[[i]][which(tmp2<=0),j] <- Min2
}


V <- rowSums(ListRandom[[i]])
dvd <- ListRandom[[i]]/V
ListRandom[[i]][which(V>=1),] <- dvd[which(V>=1),]
Vup <- rowSums(ListRandom[[i]])
VM <- 1-Vup
VM[which(VM<=0)] <- 0

V2 <- rowSums(List[[i]])
dvd2 <- List[[i]]/V2
List[[i]][which(V2>=1),] <- dvd2[which(V2>=1),]
V2up <- rowSums(List[[i]])
V2M <- 1-V2up
V2M[which(V2M<=0)] <- 0

ListRandom3D[[i]] <- cbind(ListRandom[[i]],VM)
List3D[[i]] <- cbind(List[[i]],V2M)
}



ParList <- sapply(1:length(List),function(i)DirichletempH0(List3D[[i]],ListRandom3D[[i]],
UpdatePi=UpdatePi,NumPat=NumPat,StateNames=StateNames,method=method, iter=maxround),
simplify=FALSE)

Out <- ParList

if(NumCond==2)Out <- ParList[[1]]

Out
}
