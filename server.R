library(shiny)
library(shinyFiles)
library(gdata)
library(SCPattern)
library(EBSeq)



# Define server logic for slider examples
shinyServer(function(input, output, session) {
	volumes <- c('home'="~")
    	shinyDirChoose(input, 'Outdir', roots=volumes, session=session, restrictions=system.file(package='base'))
   	 output$Dir <- renderPrint({parseDirPath(volumes, input$Outdir)})


	In <- reactive({
	print(input$Outdir)
  	outdir <- paste0("~/",input$Outdir[[1]][[2]],"/")
   	print(outdir)
	
	the.file <- input$filename$name
	if(is.null(the.file))stop("Please upload data")
	Sep=strsplit(the.file,split="\\.")[[1]]
	 if(Sep[length(Sep)]%in%c("xls"))a1=read.xls(input$filename$datapath,stringsAsFactors=F,header=TRUE, row.names=1)
	 if(Sep[length(Sep)]=="csv")a1=read.csv(input$filename$datapath,stringsAsFactors=F,header=TRUE, row.names=1)
  	if(Sep[length(Sep)]%in%c("txt","tab"))a1=read.table(input$filename$datapath,stringsAsFactors=F,header=TRUE, row.names=1)
	 Data=data.matrix(a1)

	Group.file <- input$ConditionVector$name
	if(is.null(Group.file))stop("Please upload condition file")
  	#GroupB=FALSE
  	
	Group.Sep=strsplit(Group.file,split="\\.")[[1]]
	if(Group.Sep[length(Group.Sep)]%in%c("xls"))
	GroupVIn=read.xls(input$ConditionVector$datapath,stringsAsFactors=F,header=F)
	if(Group.Sep[length(Group.Sep)]=="csv")
	GroupVIn=read.csv(input$ConditionVector$datapath,stringsAsFactors=F,header=F)
	if(Group.Sep[length(Group.Sep)]%in%c("txt","tab"))
	GroupVIn=read.table(input$ConditionVector$datapath,stringsAsFactors=F,header=F, sep="\t")
	GroupV=GroupVIn[[1]]

	if(length(GroupV)!=ncol(Data)) stop("length of the condition vector is not the same as number of cells!")

		# Compose data frame
		#input$filename$name
		List <- list(
		Input=the.file,
		GroupFile=Group.file,
		iters=input$iters, Circular=ifelse(input$Circular=="1", TRUE, FALSE),
		Cond=factor(GroupV, levels=unique(GroupV)),# follow the order they appeared
		test=input$test_buttons,
		test_details=switch(input$test_buttons,"1"="DE vs. EE (undirectional)", 
							"2"="Up vs. Down vs. EE", "3"="Up vs. Down vs. EE vs. Both direction",
							"4"="Up vs. Down"),
		RMTF=ifelse(input$RM_buttons=="1",TRUE,FALSE), 
		Dropupper=input$Dropupper,
		#PPcut=input$PPcut,
		LODNum=input$LOD, 
		Dir=outdir, 
		exExpF = paste0(outdir,input$exNormFileName,".csv"),
		exOEF = paste0(outdir,input$exListFileName,".csv"),		
		exPVF = paste0(outdir,input$exPVFileName,".csv"),
		exPVFraw = paste0(outdir,input$exPVFileName),
		PlotTF = ifelse(input$Plot_buttons=="1",TRUE,FALSE), 
		whetherLog = ifelse(input$log_whether=="1",FALSE,TRUE),
		PlotType = input$Plot_type,
		PlotF = paste0(outdir,input$exPlotFileName,".pdf"),
		PlotN = input$PlotNum
)
		# normalization and LOD
		
		if(List$test=="1"){
		Directional <- FALSE
		NumPat <- 2
		}
		else{
		  Directional <- TRUE
		  NumPat <- switch(List$test, "2"=3, "3"=4, "4"=2)  
		}
		

	 Sizes <- MedianNorm(Data)
	if(is.na(Sizes[1])){
		Sizes <- MedianNorm(Data, alternative=TRUE)
		message("alternative normalization method is applied")
	}
	 DataUse0 <- GetNormalizedMat(Data,Sizes)
 	 	
  	DataUse=DataUse0[which(apply(DataUse0,1,max)>List$LODNum),]
		# main function
  			
  	Res <- SCPTest(DataUse, List$Cond, Sizes, Circular=List$Circular, maxround=List$iters,
				 Dropout.remove=List$RMTF, Dropout.upper=List$Dropupper, LOD=List$LOD,
				 NumPat=NumPat, Directional=Directional)	
  	print("writting output...")
		browser()
	
		levs <- levels(List$Cond)
		nlevs <- nlevels(List$Cond)
		if(nlevs==2){
		PPmat <- Res$PP
		AllEE <- NA
		if("NC"%in%colnames(PPmat))AllEE <- "NC"
		if("EE"%in%colnames(PPmat))AllEE <- "EE"
		names(Res$MAP)=names(Res$maxPP)=rownames(Res$PP)
		NotEE <- names(Res$MAP)[which(Res$MAP!=AllEE)]
		NotEE.s <- names(sort(Res$maxPP[NotEE], decreasing=T))
		PPmat.sig <- cbind(NotEE.s, Res$maxPP[NotEE.s], Res$MAP[NotEE.s])
		colnames(PPmat.sig) <- c("gene", "PP most likely pattern", "most likely pattern")
		rownames(PPmat.sig) <- NotEE
		}
		if(nlevs > 2){
		PPmat <- Res$PP.all
		PPmat.sig0 <- Res$sortedlist
		NotEE.s <- rownames(Res$sortedlist)
		PPmat.sig <- cbind(NotEE.s, PPmat.sig0)
		colnames(PPmat.sig) <- c("gene", "PP most likely pattern", "most likely pattern")
		}

	write.csv(PPmat.sig,file=List$exOEF)
 	write.csv(PPmat,file=List$exPVF)
  	write.csv(DataUse0, file=List$exExpF)
		#for(i in 1:length(Res$Res))write.csv(Res$Res[[i]]$PP, file=paste0(List$exPVFraw,
		#		 "_",levs[i],"_",levs[i+1],".csv"))
	
		if(List$PlotTF){
		PN <- NULL
		if(List$PlotN!="")PN <- as.numeric(List$PlotN)
		else PN <- min(100,length(NotEE.s))
		
		if(length(PN)>0){
		pdf(List$PlotF, height=15,width=15)
		par(mfrow=c(5,4))
		for(i in 1:PN){
		if(List$PlotType%in%c("1","3")){
			if(List$whetherLog==TRUE)VioFun(NotEE.s[i],log2(DataUse0+1), 
						List$Cond, Dropout.remove=FALSE, ylab="log2(expression+1)",
						main_pre="Dropouts are shown")
			if(List$whetherLog==FALSE)VioFun(NotEE.s[i],DataUse0, 
						List$Cond, Dropout.remove=FALSE, ylab="expression",
						main_pre="Dropouts are shown")
		}

		if(List$PlotType%in%c("2","3")){
			if(List$whetherLog==TRUE)VioFun(NotEE.s[i],log2(DataUse0+1), 
				List$Cond, Dropout.remove=TRUE, ylab="log2(expression+1)",
			      	Dropout.upper=List$Dropupper, main_pre="Dropouts are not shown")
			if(List$whetherLog==FALSE)VioFun(NotEE.s[i],DataUse0, 
				List$Cond, Dropout.remove=TRUE, ylab="expression",
				Dropout.upper=List$Dropupper, main_pre="Dropouts are not shown")
		}
		}
		dev.off()
		}}

	 	List=c(List, list(Sig=PPmat.sig))	
}) 

  Act <- eventReactive(input$Submit,{
		      In()})
	# Show the values using an HTML table
  output$print0 <- renderText({
				tmp <- Act()
				str(tmp)
				paste("output directory:", tmp$Dir)
  })

	output$tab <- renderDataTable({
		tmp <- Act()$Sig
		t1 <- tmp
		print("done")
		t1
		},options = list(lengthManu = c(4,4), pageLength = 20))

#	output$done <- renderText({"Done"})
})
