#' @title Perform directional or undirectional tests across two or more conditions
#' @usage SCPTest <- function(Data, Conditions, sizeFactors,Directional = TRUE, 
#' maxround=5, NumPat=3, Circular=FALSE,Dropout.remove=FALSE, Dropout.upper=0,
#' StateNames=c("EE","Down","Up"),UpdatePi=FALSE, Print=TRUE,SmallNum=1,  Seed=10, LOD=50)
#' @param Data A data matrix contains expression values for each transcript
#'    (gene or isoform level). In which rows should be transcripts
#'    and columns should be samples. The matrix is expected to be unnormalized (raw data).
#' @param Conditions A factor indicates the condition which each sample belongs to.
#' @param sizeFactors The normalization factors. It should be a vector with lane
#'    specific numbers (the length of the vector should be the same
#'    as the number of samples, with the same order as the columns of Data).
#' @param Directional whether perform directional tests; default is TRUE.
#' If it is false, NumPat will be ignored.
#' @param Circular If Circular=TRUE, the function will also test last condition vs. first condition. Default is FALSE.
#' @param Droupout.remove,Dropout.upper For a given gene, whether cells with dropout will be 
#' ignored in the test. Defaulte is FALSE. If it is set to TRUE, values that are less or equal to
#' Dropout.upper will be ignored during the test.
#' @param maxround Number of iterations. The default value is 5.
#' @param NumPat number of patterns to consider. Could be 2, 3 or 4. Default is 3.
#' When NumPat = 3, three states (no change (EE), down and up) are considered for each 
#' transition. When NumPat = 2, two stages (down and up) are considered.
#' When NumPat = 4, four stages (no change, down, up and both directions) are considered.
#' NumPat will be ignored if Directional is FALSE. In this case the possible states are
#' no change and differentially expressed (DE).
#' @param UpdatePi default is FALSE
#' @param Print whether print the messages
#' @param SmallNum When calculating log2 expression, a small number was added to each of
#' the values. Default is 1.
#' @param Seed seed
#' @param LOD Genes with max value < LOD will be eliminated in the analysis.
#' @param method one of DirReg, VGAM and optim. It defines the method to use to estimate hyper-parameters
#' @examples mat <- matrix(exp(rnorm(100,10,1)),ncol=10)
#' Sizes <- MedianNorm(mat)
#' tmp <- SCPTest(mat, rep(c("c1","c2"),each=5), Sizes, LOD=0)
#' @author Ning Leng
#' @return Output: For each comparison, output posterior probabilities of being each state.


SCPTest <- function(
    Data, Conditions, sizeFactors,
    Directional = TRUE, Circular=FALSE,
		Dropout.remove=FALSE, Dropout.upper=0,
		maxround=5,
    NumPat=3,
    UpdatePi=FALSE, Print=TRUE,
    SmallNum=1, Seed=10, LOD=10, method="DirReg"
  ){
expect_true(NumPat%in%c(2:4))
expect_equal(ncol(Data),length(Conditions))
expect_equal(ncol(Data),length(sizeFactors))
expect_is(sizeFactors, c("numeric","integer"))
if(!is.factor(Conditions))Conditions <- factor(Conditions, levels=unique(Conditions))

if(NumPat==3) {
	if(Directional==TRUE)StateNames <- c("EE", "Down", "Up")
	if(Directional==FALSE){
		StateNames <- c("EE","DE")
		message("Directional = FALSE, only DE and EE are considered")
	}
}
if(NumPat==2){
 	if(Directional==TRUE)	StateNames <- c("Down","Up")
	if(Directional==FALSE)StateNames <- c("EE","DE")
}
if(NumPat==4) StateNames <- c("EE","Down","Up","Both")

if(Circular==TRUE){
	lv <- levels(Conditions)
	N <- length(lv)
	Which <- which(Conditions==lv[1])
	Data.which <- Data[, Which]
	colnames(Data.which) <- paste0(colnames(Data.which),"_r")
  Data <- cbind(Data, Data.which)
	lv.new <- paste0(lv[1], "_r")
  Conditions.chr <- as.character(Conditions)
	Conditions.chr2 <- c(Conditions.chr, paste0(Conditions.chr[Which],"_r"))
  Conditions <- factor(Conditions.chr2, levels=c(lv, lv.new))
	sizeFactors <- 	c(sizeFactors, sizeFactors[Which])
}
if(Directional){
  Res <- suppressWarnings(SCPairwiseDir(Data=Data, Conditions=Conditions, sizeFactors=sizeFactors,
    maxround=maxround, NumPat=NumPat, StateNames=StateNames, UpdatePi=UpdatePi,
    Print=Print, SmallNum=SmallNum, Seed=Seed, LOD=LOD, method=method,
		Dropout.remove=Dropout.remove, Dropout.upper=Dropout.upper))
}

if(!Directional){
  Res <- suppressWarnings(SCPairwiseNoDir(Data=Data, Conditions=Conditions, sizeFactors=sizeFactors,
    maxround=maxround,  StateNames=c("EE","DE"), UpdatePi=UpdatePi,
    Print=Print, SmallNum=SmallNum, Seed=Seed, LOD=LOD,
		Dropout.remove=Dropout.remove, Dropout.upper=Dropout.upper))
}

if(length(levels(Conditions))==2)Out <-  Res
else{
MAPath.list <- sapply(Res,function(i)i$MAP)
maxPP.list <- sapply(Res,function(i)i$maxPP)
if("EE"%in%colnames(Res[[1]]$PP))PPEE.list <- sapply(Res,function(i)i$PP[,"EE"])
maxPP <- apply(maxPP.list,1,prod)
EEPP <- NA
if("EE"%in%colnames(Res[[1]]$PP))EEPP <- apply(PPEE.list,1,prod)
MAPath.char <- apply(MAPath.list,1,function(i)paste0(i,collapse="-"))
names(MAPath.char) = names(maxPP) = rownames(Res[[1]]$PP)
maxPP.sort <- sort(maxPP, decreasing=TRUE)
mat.sort <- cbind(round(maxPP.sort,3), MAPath.char[names(maxPP.sort)])
colnames(mat.sort) <- c( "PP_marginal","Path")
notofi <- paste0(rep("EE",length(Res)),collapse="-")
mat.sort.use <- mat.sort[which(mat.sort[,2]!=notofi),]

# marginal for all
all.pathname <- StateNames
for(i in 1:(length(Res)-1))
		 all.pathname<- as.vector(outer(all.pathname,StateNames,function(i,j)paste(i,j,sep="-")))
path.mat <- t(sapply(1:length(all.pathname),function(i)strsplit(all.pathname[i],split="-")[[1]]))
path.all <- sapply(1:length(all.pathname),function(i)
						apply(sapply(1:length(Res),function(j)Res[[j]]$PP[,path.mat[i,j]]),1,prod))
rownames(path.mat)=colnames(path.all)=all.pathname
margin.maxval <- apply(path.all,1,max) 
margin.max <- all.pathname[apply(path.all,1,which.max)]
names(margin.max) <- names(margin.maxval)
Out <- list(sortedlist=mat.sort.use, Path=MAPath.char, Path.mat=MAPath.list, maxPP=maxPP,
						maxmarginPP=margin.maxval, maxmarginPath=margin.max,
						EEPP=EEPP,PP.all=path.all, path.mat=path.mat,all.pathname = all.pathname, Res=Res)
}

Out
}
