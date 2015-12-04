# SCPattern
Identify and characterize expression changes in a single cell RNA-seq experiment with ordered conditions

## Installation
The SCPattern GUI requires the following packages : shiny, shinyFiles, SCPattern

To install the shiny packages, in R run:

> install.packages('shiny')

> install.packages('shinyFiles')

SCPattern R package and its vignette could be found at https://github.com/lengning/SCPattern/tree/master/package

To install these two packages, in bash run 

> R CMD INSTALL SCPattern_0.0.3.tar.gz

Or install locally via R GUI / Rstudio

## Run the app

In R, run:

> library(shiny)

> runGitHub('lengning/SCPattern')

![Screenshot](https://github.com/lengning/SCPattern/blob/master/figs/SCPattern_screenshot.png)


## Input files

The first input file should be the expression matrix. 
Rows are the genes and columns are the cells.
Currently the program only takes csv files or tab delimited file.
The input file will be treated as a tab delimited file if the suffix is not '.csv'.


The second input file is the condition vector. It could be csv or tab delimited file. The file should contain
1 column. The i th component represents the condition that cell i belongs to. The length of the condition vector
should be the same as the number of cells in the first input file.
Two or more conditions are expected. 

### Two conditions

### Multiple conditions

An example input file **example_sc_data** could be found at https://github.com/lengning/SCPattern/tree/master/example_data/   
- This file contains 1000 genes and 300 cells


An example condition vector file **example_sc_conditionvector** could be found at https://github.com/lengning/SCPattern/tree/master/example_data/
- This file contains 300 elements; the 300 cells are from 5 different conditions- 

Another set of expression input file and condition vector file **example_sc_data_twocond** and **example_sc_conditionvector_twocond** can also be found under https://github.com/lengning/SCPattern/tree/master/example_data/
- In this set of files, there are 120 cells from two conditions.

## Note

The 'create new folder' button in the output folder selection pop-up is disfunctional right now
