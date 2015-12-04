# SCPattern
Identify and characterize expression changes in a single cell RNA-seq experiment with ordered conditions

## 1. Installation
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


## 2. Input files

The first input file should be the expression matrix. 
Rows are the genes and columns are the cells.
Currently the program only takes csv files or tab delimited file.
The input file will be treated as a tab delimited file if the suffix is not '.csv'.


The second input file is the condition vector. It could be csv or tab delimited file. The file should contain
1 column. The i th component represents the condition that cell i belongs to. The length of the condition vector
should be the same as the number of cells in the first input file.
Two or more conditions are expected. 

### 2.1 Two conditions
An example input file **example_sc_data_twocond** could be found at https://github.com/lengning/SCPattern/tree/master/example_data/   
- This file contains 1000 genes and 120 cells

An example condition vector file **example_sc_conditionvector_twocond** could be found at https://github.com/lengning/SCPattern/tree/master/example_data/
- This file contains 120 elements; the 120 cells are from 2 different conditions- 

### 2.2 Multiple conditions
An example input file **example_sc_data** could be found at https://github.com/lengning/SCPattern/tree/master/example_data/   
- This file contains 1000 genes and 300 cells

An example condition vector file **example_sc_conditionvector** could be found at https://github.com/lengning/SCPattern/tree/master/example_data/
- This file contains 300 elements; the 300 cells are from 5 different conditions- 

## 3. Features

number of iteration
## Note

The 'create new folder' button in the output folder selection pop-up is disfunctional right now
