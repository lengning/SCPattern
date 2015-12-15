#' @title Calculate non-directional KS-statistics of pair-wise comparisons
#' @usage KSNoDir(SCMat, Cond, pval=FALSE,Dropout.remove=FALSE, Dropout.upper=0)
#' @param SCMat A data matrix contains expression values for each transcript
#'    (gene or isoform level). In which rows should be transcripts
#'    and columns should be samples. The matrix is expected to be normalized
#'    prior to the analysis.
#' @param Cond A factor indicates the condition which each sample belongs to.
#' @param pval whether output KS test raw p values
#' @param Droupout.remove,Dropout.upper For a given gene, whether cells with dropout will be 
#' ignored in the test. Defaulte is FALSE. If it is set to TRUE, values that are less or equal to
#' Dropout.upper will be ignored during the test.
#' @return Output contains KS statistics (non-directional) 
#' for each pairwise comparison. 
#' Note the order of conditions are defined by the order in the Cond factor.
#' @author Ning Leng
#' @examples mat <- matrix(rnorm(40), ncol=10)
#' tmp1 <- KSNoDir(mat, as.factor(rep(c("c1","c2"),each=5))


KSNoDir <- function(SCMat, Cond, pval=FALSE, Dropout.remove=FALSE, Dropout.upper=0){
SCMatUselog <- SCMat
SCCondUse <- Cond
if(Dropout.remove==FALSE)
StatAllgene <- sapply(1:nrow(SCMatUselog),function(k){
	tmpgene <- rownames(SCMatUselog)[k]
	statout <- sapply(1:(nlevels(SCCondUse)-1),function(i)
	ks.test(SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i])],
				SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i+1])]
				)$statistic)
statout})
  if(Dropout.remove==TRUE)
  StatAllgeneGt <- sapply(1:nrow(SCMatUselog),function(k){
  tmpgene <- rownames(SCMatUselog)[k]
  statout <- sapply(1:(nlevels(SCCondUse)-1),function(i){
  t1 <- SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i])]
  t2 <- SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i+1])]
  t3 <- t1[which(t1>Dropout.upper)]
  t4 <- t2[which(t2>Dropout.upper)]
  if(length(t3)>3 & length(t4) >3)out <- ks.test(t3,t4)$statistic
  else out <- 0})
  statout})



if(nlevels(SCCondUse)>2)StatMat <- t(StatAllgene)
if(nlevels(SCCondUse)==2)StatMat <- matrix(StatAllgene,nrow=nrow(SCMatUselog))

rownames(StatMat) <- rownames(SCMatUselog)
StatMatMean <- rowMeans(StatMat)
StatMatMeanSort <- sort(StatMatMean,decreasing=TRUE)
StatMatMeanRank <- length(StatMatMean)+1-rank(StatMatMean)

StatMatSort <- sapply(1:(nlevels(SCCondUse)-1),function(j){
tmp <- abs(StatMat[,j])
t2 <- names(sort(tmp,decreasing=TRUE))
})

Out <- list(Test = StatMat, TestSort=StatMatSort)


if(pval==TRUE){
if(Dropout.remove==FALSE)
StatAllgene <- sapply(1:nrow(SCMatUselog),function(k){
	tmpgene <- rownames(SCMatUselog)[k]
	statout <- sapply(1:(nlevels(SCCondUse)-1),function(i)
	ks.test(SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i])],
				SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i+1])]
				)$p.value)
statout})
  if(Dropout.remove==TRUE)
  StatAllgeneGt <- sapply(1:nrow(SCMatUselog),function(k){
  tmpgene <- rownames(SCMatUselog)[k]
  statout <- sapply(1:(nlevels(SCCondUse)-1),function(i){
  t1 <- SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i])]
  t2 <- SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i+1])]
  t3 <- t1[which(t1>Dropout.upper)]
  t4 <- t2[which(t2>Dropout.upper)]
  if(length(t3)>3 & length(t4) >3)out <- ks.test(t3,t4)$p.value
  else out <- 0})
  statout})
if(nlevels(SCCondUse)>2)StatMat <- t(StatAllgene)
if(nlevels(SCCondUse)==2)StatMat <- matrix(StatAllgene,nrow=nrow(SCMatUselog))
rownames(StatMat) <- rownames(SCMatUselog)
StatMatMean <- rowMeans(StatMat)
StatMatMeanSort <- sort(StatMatMean,decreasing=TRUE)
StatMatMeanRank <- length(StatMatMean)+1-rank(StatMatMean)

StatMatSort <- sapply(1:(nlevels(SCCondUse)-1),function(j){
tmp <- abs(StatMat[,j])
t2 <- names(sort(tmp))
})

pOut <- list(
pvals=StatMat,  pvalsSort=StatMatSort)

Out <- c(Out,pOut)
}

Out
}
