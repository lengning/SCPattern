\name{SCPatternExData}
\alias{SCPatternExData}
\docType{data}
\title{
	Simulated gene level data set with 5 ordered conditions
}
\description{
	'SCPatternExData' gives the gene level simulated data with 5 ordered conditions, 60 cells for each condition.
	The data set was simulated following the Negative Binomial distribution. The parameters of each gene (mean and overdispersion)
	were sampled from the empirical estimates from an empirical RNA-Seq data set from Thomson lab at Morgridge Institute for 
	Research.  
}
\format{
	SCPatternExData is a matrix with 1000 genes (rows) and 300 cells (columns). 
}
\seealso{
}
\examples{
data(SCPatternExData)
str(SCPatternExData)
}
\keyword{datasets}




