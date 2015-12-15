#' @title Perform non-directional tests across two or more conditions
#' @usage SCPairwiseNoDir <- function(Data, Conditions, sizeFactors,
#' Dropout.remove=FALSE, Dropout.upper=0,maxround=5,
#' StateNames=c("EE","DE"),UpdatePi=TRUE, Print=TRUE,SmallNum=1,  Seed=10, LOD=50)
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
#' @param StateNames names of the states
#' @param UpdatePi default is TRUE
#' @param Print whether print the messages
#' @param SmallNum When calculating log2 expression, a small number was added to each of
#' the values. Default is 1.
#' @param Seed seed
#' @param LOD Genes with max value < LOD will be eliminated in the analysis.
#' @examples mat <- matrix(exp(rnorm(100,10,1)),ncol=10)
#' Sizes <- MedianNorm(mat)
#' tmp <- SCPairwiseNoDir(mat, rep(c("c1","c2"),each=5), Sizes, LOD=0)
#' @author Ning Leng
#' @return Output: For each comparison, output posterior probabilities of being each state.


SCPairwiseNoDir <- function(
		Data, Conditions, sizeFactors,
		Dropout.remove=FALSE, Dropout.upper=0,maxround=5,
		StateNames=c("EE","DE"),
		UpdatePi=TRUE, Print=TRUE, 
		SmallNum=1, Seed=10, LOD=50
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

DataList.unlist.dvd.log <- log2(DataList.unlist.dvd + SmallNum)

###################################
# Calculate KS stat for all comparisons
###################################
TestRes <- KSNoDir(DataList.unlist.dvd.log, Conditions,Dropout.remove=Dropout.remove, 
									 Dropout.upper=log2(Dropout.upper+SmallNum))
List <- TestRes$Test


# shuffle
Which <- sapply(1:NumCond,function(i)which(Conditions==CondLevels[i]),simplify=FALSE)
Shuffle <- sapply(Which,function(i)sample(i,length(i)))
ShuffleV <- unlist(Shuffle)
DataShuffle <- DataList.unlist.dvd.log[,ShuffleV]

ListRandom <- KSNoDir(DataShuffle, Conditions,Dropout.remove=Dropout.remove, Dropout.upper=Dropout.upper)$Test


#################################
# impute 0s

List3D=ListRandom3D=vector("list",ncol(List))
for(i in 1:ncol(List)){
t0 <- List[,i]
t0[which(t0<0)] <- 0
t0[which(t0>1)] <- 1
List3D[[i]] <- t0
names(List3D[[i]]) <- rownames(DataList.unlist.dvd.log)
t1 <- ListRandom[,i]
t1[which(t1<0)] <- 0
t1[which(t1>1)] <- 1
ListRandom3D[[i]] <- t1
names(ListRandom3D[[i]]) <- rownames(DataList.unlist.dvd.log)
}



ParList <- sapply(1:length(List3D),function(i)BetaempH0(List3D[[i]],ListRandom3D[[i]],
UpdatePi=UpdatePi,StateNames=StateNames),
simplify=FALSE)

Out <- ParList

if(NumCond==2)Out <- ParList[[1]]

Out
}
