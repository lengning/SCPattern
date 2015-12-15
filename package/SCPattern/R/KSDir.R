#' @title Calculate directional KS-statistics of pair-wise comparisons
#' @usage KSDir(SCMat, Cond, pval=FALSE, Dropout.remove=FALSE, Dropout.upper=0)
#' @param SCMat A data matrix contains expression values for each transcript
#'    (gene or isoform level). In which rows should be transcripts
#'    and columns should be samples. The matrix is expected to be normalized
#'    prior to the analysis.
#' @param Cond A factor indicates the condition which each sample belongs to.
#' @param pval whether output KS test raw p values
#' @param Droupout.remove,Dropout.upper For a given gene, whether cells with dropout will be 
#' ignored in the test. Defaulte is FALSE. If it is set to TRUE, values that are less or equal to
#' Dropout.upper will be ignored during the test.
#' @return Output contains KS statistics (up-regulate and down-regulate) 
#' for each pairwise comparison. 
#' Note the order of conditions are defined by the order in the Cond factor.
#' @author Ning Leng
#' @examples mat <- matrix(rnorm(40), ncol=10)
#' tmp1 <- KSDir(mat, as.factor(rep(c("c1","c2"),each=5))


KSDir<-function(SCMat, Cond, pval=FALSE, Dropout.remove=FALSE, Dropout.upper=0){
  if(ncol(SCMat)!=length(Cond))stop("Error: length of condition vector is different from number of columns of the input data")
  SCMatUselog <- SCMat
  SCCondUse <- Cond
  if(Dropout.remove==FALSE)
  StatAllgeneGt <- sapply(1:nrow(SCMatUselog),function(k){
	tmpgene <- rownames(SCMatUselog)[k]
	statout <- sapply(1:(nlevels(SCCondUse)-1),function(i)
	ks.test(SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i])],
				SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i+1])],
				alternative="greater")$statistic)
  statout})
  if(Dropout.remove==TRUE)
  StatAllgeneGt <- sapply(1:nrow(SCMatUselog),function(k){
  tmpgene <- rownames(SCMatUselog)[k]
  statout <- sapply(1:(nlevels(SCCondUse)-1),function(i){
  t1 <- SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i])]
  t2 <- SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i+1])]
  t3 <- t1[which(t1>Dropout.upper)]
  t4 <- t2[which(t2>Dropout.upper)]
  if(length(t3)>3 & length(t4) >3)out <- ks.test(t3,t4,alternative="greater")$statistic
  else out <- 0})
  statout})


  if(Dropout.remove==FALSE)
  StatAllgeneLess <- sapply(1:nrow(SCMatUselog),function(k){
  tmpgene <- rownames(SCMatUselog)[k]
  statout <- 
  sapply(1:(nlevels(SCCondUse)-1),function(i)
	ks.test(SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i])],
					SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i+1])],
					alternative="less")$statistic)
  statout})
  if(Dropout.remove==TRUE)
  StatAllgeneLess <- sapply(1:nrow(SCMatUselog),function(k){
  tmpgene <- rownames(SCMatUselog)[k]
  statout <- sapply(1:(nlevels(SCCondUse)-1),function(i){
  t1 <- SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i])]
  t2 <- SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i+1])]
  t3 <- t1[which(t1>Dropout.upper)]
  t4 <- t2[which(t2>Dropout.upper)]
  if(length(t3)>3 & length(t4) >3)out <- ks.test(t3,t4,alternative="less")$statistic
  else out <- 0})
  statout})



if(nlevels(SCCondUse)>2)StatMatGt <- t(StatAllgeneGt)
if(nlevels(SCCondUse)==2)StatMatGt <- matrix(StatAllgeneGt,nrow=nrow(SCMatUselog))

rownames(StatMatGt) <- rownames(SCMatUselog)
StatMatGtMean <- rowMeans(StatMatGt)
StatMatGtMeanSort <- sort(StatMatGtMean,decreasing=TRUE)
StatMatGtMeanRank <- length(StatMatGtMean)+1-rank(StatMatGtMean)

StatMatGtSort <- sapply(1:(nlevels(SCCondUse)-1),function(j){
tmp <- abs(StatMatGt[,j])
t2 <- names(sort(tmp,decreasing=TRUE))
})

if(nlevels(SCCondUse)>2)StatMatLess <- t(StatAllgeneLess)
if(nlevels(SCCondUse)==2)StatMatLess <- matrix(StatAllgeneLess,nrow=nrow(SCMatUselog))
rownames(StatMatLess) <- rownames(SCMatUselog)
StatMatLessMean <- rowMeans(StatMatLess)
StatMatLessMeanSort <- sort(StatMatLessMean,decreasing=TRUE)
StatMatLessMeanRank <- length(StatMatLessMean)+1-rank(StatMatLessMean)

StatMatLessSort <- sapply(1:(nlevels(SCCondUse)-1),function(j){
	 tmp <- abs(StatMatLess[,j])
	 t2 <- names(sort(tmp,decreasing=TRUE))
	 })

CB <- sapply(1:(nlevels(SCCondUse)-1),function(j){
tmp <- cbind(StatMatGt[,j],StatMatLess[,j])
tmp},simplify=FALSE)

Out <- list(List=CB,
Greater=StatMatGt, Less=StatMatLess, GreaterSort=StatMatGtSort, LessSort=StatMatLessSort)


if(pval==TRUE){
  if(Dropout.remove==FALSE)
  StatAllgeneGt <- sapply(1:nrow(SCMatUselog),function(k){
	tmpgene <- rownames(SCMatUselog)[k]
	statout <- sapply(1:(nlevels(SCCondUse)-1),function(i)
	ks.test(SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i])],
				SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i+1])],
				alternative="greater")$p.value)
  statout})
  if(Dropout.remove==TRUE)
  StatAllgeneGt <- sapply(1:nrow(SCMatUselog),function(k){
  tmpgene <- rownames(SCMatUselog)[k]
  statout <- sapply(1:(nlevels(SCCondUse)-1),function(i){
  t1 <- SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i])]
  t2 <- SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i+1])]
  t3 <- t1[which(t1>Dropout.upper)]
  t4 <- t2[which(t2>Dropout.upper)]
  if(length(t3)>3 & length(t4) >3)out <- ks.test(t3,t4,alternative="greater")$p.value
  else out <- 1})
  statout})

  if(Dropout.remove==FALSE)
  StatAllgeneLess <- sapply(1:nrow(SCMatUselog),function(k){
  tmpgene <- rownames(SCMatUselog)[k]
  statout <- sapply(1:(nlevels(SCCondUse)-1),function(i)
	ks.test(SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i])],
					SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i+1])],
					alternative="less")$p.value)
statout})
  if(Dropout.remove==TRUE)
  StatAllgeneLess <- sapply(1:nrow(SCMatUselog),function(k){
  tmpgene <- rownames(SCMatUselog)[k]
  statout <- sapply(1:(nlevels(SCCondUse)-1),function(i){
  t1 <- SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i])]
  t2 <- SCMatUselog[tmpgene,which(SCCondUse==levels(SCCondUse)[i+1])]
  t3 <- t1[which(t1>Dropout.upper)]
  t4 <- t2[which(t2>Dropout.upper)]
  if(length(t3)>3 & length(t4) >3)out <- ks.test(t3,t4,alternative="less")$p.value
  else out <- 1})
statout})

if(nlevels(SCCondUse)>2)StatMatGt <- t(StatAllgeneGt)
if(nlevels(SCCondUse)==2)StatMatGt <- matrix(StatAllgeneGt,nrow=nrow(SCMatUselog))
rownames(StatMatGt) <- rownames(SCMatUselog)
StatMatGtMean <- rowMeans(StatMatGt)
StatMatGtMeanSort <- sort(StatMatGtMean,decreasing=TRUE)
StatMatGtMeanRank <- length(StatMatGtMean)+1-rank(StatMatGtMean)

StatMatGtSort <- sapply(1:(nlevels(SCCondUse)-1),function(j){
tmp <- abs(StatMatGt[,j])
t2 <- names(sort(tmp))
})

if(nlevels(SCCondUse)>2)StatMatLess <- t(StatAllgeneLess)
if(nlevels(SCCondUse)==2)StatMatLess <- matrix(StatAllgeneLess,nrow=nrow(SCMatUselog))
rownames(StatMatLess) <- rownames(SCMatUselog)
StatMatLessMean <- rowMeans(StatMatLess)
StatMatLessMeanSort <- sort(StatMatLessMean,decreasing=TRUE)
StatMatLessMeanRank <- length(StatMatLessMean)+1-rank(StatMatLessMean)

StatMatLessSort <- sapply(1:(nlevels(SCCondUse)-1),function(j){
	 tmp <- abs(StatMatLess[,j])
	 t2 <- names(sort(tmp))
	 })

CB <- sapply(1:(nlevels(SCCondUse)-1),function(j)cbind(StatMatGt[,j],StatMatLess[,j]),simplify=FALSE)
pOut <- list(p.List=CB,
p.Greater=StatMatGt, p.Less=StatMatLess, p.GreaterSort=StatMatGtSort, p.LessSort=StatMatLessSort)

Out <- c(Out,pOut)
}

Out
}
