# SCPattern
Identify and characterize expression changes in a single cell RNA-seq experiment with ordered conditions

## Installation
The SCPattern GUI requires the following packages : shiny, shinyFiles, SCPattern

To install the shiny packages, in R run:

> install.packages("shiny")

> install.packages("shinyFiles")

SCPattern R package and its vignette could be found at https://github.com/lengning/SCPattern/tree/master/package

To install these two packages, in bash run 

> R CMD INSTALL SCPattern_0.0.2.tar.gz


## Run the app

In R, run:

> library(shiny)

> runGitHub("lengning/SCPattern")
