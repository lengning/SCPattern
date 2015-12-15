#' @title Plot violinplot
#' @usage VioFun(gene, Data, Conditions, Col, main_pre=NULL, Dropout.remove=FALSE, Dropout.upper=0, namecex=1, ylab="Expression", Log="")
#' @param gene gene name of interest (character)
#' @param Data input data, gene-by-sample matrix. Normalized expression is expected
#' @param Conditions a vector of characters that indicates which condition the cells belong to
#' @param Col colors
#' @param main_pre prefix for main
#' @param namecex cex for x axis names 
#' @param ylab y axis label
#' @param Log whether plot in log scale
#' @param Droupout.remove,Dropout.upper For a given gene, whether cells with dropout will be 
#' ignored in the test. Defaulte is FALSE. If it is set to TRUE, values that are less or equal to
#' Dropout.upper will be ignored during the test.
#' @return violin plot 
#' @examples Data<- abs(matrix(rnorm(20,100,2), nrow=2))
#' rownames(Data) <- paste0("gene",1:2)
#' cond <- factor(rep(c("a","b"),each=5))
#' VioFun("gene1",Data, cond)
#' @author Ning Leng

VioFun=function(gene, Data, Conditions, Col=NULL, main_pre=NULL, 
			Dropout.remove=FALSE, Dropout.upper=0, namecex=1,
			ylab="Expression", Log=""){
	SCMat <- Data
	SCCond <- Conditions
	if(!is.factor(SCCond))SCCond <- factor(SCCond, levels=unique(SCCond))
	expect_equal(ncol(SCMat), length(SCCond))
	mainpaste <- main_pre
  if(is.null(Col))Col <- rainbow(nlevels(SCCond))
	tmpdata <- SCMat[gene,]
	tmpcond <- SCCond
	if(Dropout.remove) {
		wch <- which(tmpdata>Dropout.upper)
		tmpdata <- tmpdata[wch]
		tmpcond <- SCCond[wch]
	}
  SCDataVio <- sapply(1:nlevels(SCCond),function(k)
        tmpdata[tmpcond==levels(SCCond)[k]],simplify=F)
  names(SCDataVio) <- paste0(levels(SCCond)," (",sapply(SCDataVio,length),")")

  tmpmin <- 0
  if(Log=="y")tmpmin <- 1
  plot(1,1,col="white",xlim=c(0,length(SCDataVio)+1),
  ylim=c(tmpmin,(max(SCMat[gene,])+.1)*1.2),
  xaxt="none",ylab=ylab, xlab="",main=paste(gene, mainpaste),
  log=Log)
  mtext(side=1,at=1:length(SCDataVio),paste0(names(SCDataVio)," "),
las=3,cex=namecex)
  for(j in 1:length(SCDataVio))
  if(max(SCDataVio[[j]])>0)vioplot(SCDataVio[[j]],at=j, add=T,col=Col[j],colMed=NULL)
}



