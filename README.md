
# SCPattern
Identify and characterize expression changes in a single cell RNA-seq experiment with ordered conditions

## 1. Installation
The SCPattern GUI requires the following packages : shiny, shinyFiles, SCPattern

To install the shiny packages, in R run:

> install.packages('shiny')

> install.packages('shinyFiles')

SCPattern R package and its vignette could be found at https://github.com/lengning/SCPattern/tree/master/package

To install SCPattern, in R run: 

> install.packages("devtools")

> library(devtools)

> install_github("lengning/SCPattern/package/SCPattern")

Or install locally.

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

## 3. Customize options
- The number of iteration: Default is 5. 
-	Testing types: User can choose the test type they want to perform. For simple DE identification, one can choose "DE vs. EE"; for directional test, one can choose "Up vs. Down vs. EE" or "Up vs. Down vs. EE vs. Both direction" or "Up vs. Down". Depending on what you choose, result provides posterior probabilities of all possible combinations. For example, when we have 5 time points, option "Up vs. Down" will give you 2^4=16 possible combinations.
- Ignore dropouts: If you choose "Ignore", when performing DE analysis for a gene, any cells with zero count would not be considered in the analysis.
-	The droupout is defined as value <= a: User may set the dropout cutoff, default is 0.
- Lower limit of detection (max value): if it is set as m, genes whose max value is less than m will be ignored
- Circular: whether test the last condition vs. the first condition in addition to the pairwise tests.
- Output directory, will be set as home directory (~/) if it is empty.
- Output file name for the normalized expression.
- Output file name for the detected DE genes.
- Output file name for the posterior probability (PP) matrix.
-	Whether plot DE genes or not.
- Whether show zeros in plots.
-	Whether take a log scale in plot.
-	Number of genes to plot. If it is not specified, top 100 genes will be included in the output plots.
-	Output file name for the gene plots.

## 4. Outputs
Three (Four) files will be generated:
-	genes.csv: genes are listed with their corresponding PPs and the most likely pattern. Genes are sorted by PP.
- PPs.csv: For each gene (row), PPs are shown for all possible combination of patterns (column).
- normalized.csv: Normalized expression matrix with genes in row and cells in column. 
- Plots.pdf: This file will be generated only when the user chooses to plot top genes. In each plot, x-axis shows conditions and y-axis shows expression. Genes are sorted by their PP. 
 
## Note

The 'create new folder' button in the output folder selection pop-up is disfunctional right now

## License
This project is licensed under the terms of the Apache License 2.0

