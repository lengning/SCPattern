#' @title Estimate parameters of a dirichlet distribution
#' @param ParamPool parameters
#' @param InputPool input values
#' @return sum likelihood
#' @author Ning Leng

DirichletOptim<-function(ParamPool, InputPool){
	val=ParamPool
	In=InputPool
	al=val[1]
	be=val[2]
	par=c(al,al,be)
	Data=In[[1]]
	Num=-100
	library(gtools)
	Res=ddirichlet(Data,par)
	ResLog0=log(Res)
	ResLog=ifelse(ResLog0==-Inf, -100, ResLog0)
	-sum(ResLog)
}

