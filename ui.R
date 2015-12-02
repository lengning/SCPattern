library(shiny)
library(shinyFiles)
library(gdata)
options(shiny.maxRequestSize=500*1024^2) 
# Define UI for slider demo application
shinyUI(pageWithSidebar(

  #  Application title
  headerPanel("SCPattern"),

  # Sidebar with sliders that demonstrate various available options
  sidebarPanel(width=9,
    # file
    
		fileInput("filename", label = "File input (support .csv, .xls, .txt, .tab)"),

		# grouping vector
		fileInput("ConditionVector", label = "Condition vector \n file name (support .csv, .xls, .txt, .tab)"),
	
		column(3,


				numericInput("iters",
				label = "number of iteration",
							value = 5),

				# Normalization
				radioButtons("test_buttons",
						label = "test between",
					choices = list("DE vs. EE (undirectional)" = 1,
							"Up vs. Down vs. EE" = 2,
							"Up vs. Down vs. EE vs. Both direction"=3,
							"Up vs. Down"=4),
								selected = 2),								

				radioButtons("RM_buttons",
						label = "Ignor dropouts?",
						 choices = list("Ignore" = 1,
									 "Do not ignore" = 2),
												 selected = 2),

				numericInput("Dropupper",
				label = "The dropout is defined as values < =",
												value = 0)
				),								
		
		column(3,
				# LOD
				numericInput("LOD",
				label = "Lower limit of detection (max value)",
										value = 10),

				radioButtons("Circular",
				label = "Circular?",
				 choices = list("Yes" = 1,
							 "No" = 2),
										 selected = 2)#,				 
				# Num permutation
				#numericInput("PPcut",
				#label = "PP cutoff",
				#				value = 0.5)
		),

		column(3,
		    # output dir
	    	shinyDirButton('Outdir', 'output folder select', 'Please select a folder'),
		
				# export normalzied matrix
				textInput("exNormFileName", 
				label = "Export file name - normalized expression matrix", 
									value = "normalized"),
					
				# export gene list
				textInput("exListFileName", 
				label = "Export file name - DE gene list", 
									value = "genes"),
				# export PP
				textInput("exPVFileName", 
				label = "Export file name - PP matrix", 
									value = "PPs")
	),

	column(3,
		    radioButtons("Plot_buttons",
		    label = "Plot top genes?",
				 choices = list("Yes" = 1,
					        "No" = 2),
					         selected = 1),								
		
		    radioButtons("Plot_type",
		    label = "Show 0s in plots?",
				 choices = list("Yes" = 1,
					        "No" = 2,
						"In both way" = 3),
					         selected = 1),

        	radioButtons("log_whether",
     		label = "Plot in log scale?",
					  choices = list("No" = 1,
					       "log2(expression + 1)" = 2),
					         selected = 1),
				# num genes to plot
				textInput("PlotNum", 
				label = "Number of genes to plot (if not specified, top 100 genes will be plotted)", 
								value = ""),
				# plot name
				textInput("exPlotFileName", 
				label = "Export file name for the plots?", 
								value = "Plots")
		),

		actionButton("Submit","Submit for processing")
),

	
						

  # Show a table summarizing the values entered
  mainPanel(
    h4(textOutput("print0")),
		#tableOutput("values")
		dataTableOutput("tab")
  )
))
