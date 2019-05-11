rm(list=ls())
getwd()
setwd("C:/Users/Jose/Documents/Ramapo/Bioinformatics/FinalProject/Leader_Genes_Approach/")

## Jose Carrillo - Shiny Application - Advanced Bioinformatics - Spring 2019

# This is the front-end for the Leader Genes Approach application consisting of two main panels (or "pages") controlled as conditional panels:
#	1) The title panel/page allows the user to either enter their own search term OR import a .csv file of genes, it is hidden afterwards
#	2) The data panel/page is displayed after a loading screen, it consists of multiple tabs that hold relevant data and a final tab for exporting data

library("DT")
library("ggplot2")
library("reticulate")
library("shiny")
library("shinyjs")
library("STRINGdb", lib.loc="library")

source("scripts/processing.R")
source_python("scripts/mining.py")  # Must be in the same working directory

#---------- VARIABLE INITIALIZATION ----------#

term <- ""
genes_interactions <- ""
genes_scores <- ""
leader_genes <- ""

#---------- SHINY ----------#

ui <- fluidPage (
	useShinyjs()

	, tags$style (  # Visually appealing! Very aesthetic! I use this combo for all of my presentations!
		HTML("@import url('//fonts.googleapis.com/css?family=Raleway');")
		, HTML("@import url('//fonts.googleapis.com/css?family=Merriweather');")
		, HTML(".shiny-notification { height: 100px; width: 800px; position:fixed; top: calc(50% - 50px); left: calc(50% - 400px); font-size: 250%; text-align: center; }")
	)
	
	# This conditional panel is the "title page" and is hidden once the user input is taken
	, div( id="title_animation", conditionalPanel (
		condition = "!output.hide_title"
		, fluidRow (
			style="padding:200px;"
			, fluidRow (  # Page heading and a brief description of the tool
				style="padding:20px;"
				, column(width=12, align="center", 
					h1(style="font-family: Merriweather;", "Leader Genes Approach: An ", em("ab initio"), " algorithm for predicting functional genes of complex diseases"))
				, column(width=12, align="center", 
					h3(style="font-family: Raleway;", "The Leader Genes Approach (LGA) utilizes a novel algorithm to find a preliminary set of genes that are functionally associated with a particular disease. Data is mined from the National Center for Biotechnology Information (NCBI) and interactions are taken from the Search Tool for the Retrieval of Interacting Genes/Proteins (STRING) database. Note that the LGA is set to enforce a strict 95% confidence for interactions.")
					, h3(style="font-family: Raleway;", "Please enter a well-known disease in the search bar to begin!"))
			)

			, fluidRow (  # Text input for search term and enter button to initialize
				column(width=7, align="right", style="font-family: Raleway;",
					textInput(inputId="search_term.text", label="", placeholder="Periodontitis"))
				, column(width=5, align="left", style="font-family: Raleway; margin-top: 20px;", class="form-group shiny-input-container", 
					actionButton(inputId="search_term.button", label="Enter"))
			)

			, fluidRow (  # File input for a set of genes and import button to initialize
				column(width=7, align="right", style="font-family: Raleway;",
					fileInput(inputId="import.file", label="", multiple=FALSE, accept=c(".csv"), buttonLabel="Browse...", placeholder="Import a set of genes"))
				, column(width=5, align="left", style="font-family: Raleway; margin-top: 20px;", class="form-group shiny-input-container", 
					actionButton(inputId="import.button", label="Import"))
			)
		)
	))

	# This conditional panel contains all of the processed data and is initially hidden until the user input is taken
	, conditionalPanel (
		condition = "output.hide_title"
		, tabsetPanel (id=NULL, selected=NULL, type="tabs"  # Each tab contains a page of relevant data

			, tabPanel ("Leader Genes"  # Shows the top 10 leader genes of the STRING network
				, fluidRow (
					column(width=12, align="center", style="font-family: Merriweather; padding: 20px;"
						, h2("Leader Genes"))
					, column(width=12, align="center", style="font-family: Raleway; margin-bottom: 20px;"
						, h3("The top ten leader genes are listed below. The scores were calculated based on the strength and number of interactions of the STRING network."))
					, column(width=12, align="left", style="font-family: Raleway;"
						, DT::dataTableOutput(outputId="leader_genes_dt", width="auto", height="auto"))
				)
			)


			, tabPanel ("Network"  # Shows the actual STRING network as well as all interactions between genes
				, fluidRow (
					column(width=12, align="center", style="font-family: Merriweather; padding: 20px;"
				 		, h2("STRINGdb Gene Network and Weighted Interactions"))
					, column(width=12, align="center", style="font-family: Raleway; margin-bottom: 20px;"
						, h3("There are seven main evidence channels with which STRING derives its scores from: 
								(i) ", strong(em(span(style="color:magenta;", "Experiments"))), " - Evidence from lab experiment data (biochemical, biophysical, genetic) from the IMEx consortium and BioGRID; 
								(ii) ", strong(em(span(style="color:cyan;", "Database"))), " - Evidence asserted by a human expert curator and imported from pathway databases; 
								(iii) ", strong(em(span(style="color:yellowgreen;", "Textmining"))), " - Evidence from mentions of protein names in all PubMed abstracts and in other text collections (OMIM, SGD); 
								(iv) ", strong(em(span(style="color:black;", "Co-expression"))), " - Evidence from gene expression data that has been normalized, pruned and correlated; 
								(v) ", strong(em(span(style="color:limegreen;", "Neighborhood"))), " - Evidence from genome-based prediction, but mostly relevant for Bacteria and Archaea rather than humans; 
								(vi) ", strong(em(span(style="color:red;", "Fusion"))), " - An additional association score is given when a pair of proteins' respective orthologs have fused into a single, protein-coding gene in at least one organism; and 
								(vii) ", strong(em(span(style="color:blue;", "Co-occurrence"))), " - Evidence from the phylogenetic distribution of orthologs of all proteins in a given organism."))
					, column(width=12, align="center"
						, plotOutput(outputId="network_plot", width=1800, height=1350))
					, column(width=12, align="left", style="font-family: Raleway; padding: 20px"
						, DT::dataTableOutput(outputId="genes_interactions_dt", width="auto", height="auto"))
				)
			)
			

			, tabPanel ("Clusters" 	# Shows subnetworks of the STRING network, or the clusters, as well as the associated genes and their subnetwork global scores
				, fluidRow(
					column(width=12, align="center", style="font-family: Merriweather; padding: 20px;"
						, h2("Clusters"))
					, column(width=12, align="center", style="font-family: Raleway; margin-bottom: 20px;"
						, h3("The first five clusters of the STRING network along with their gene scores are listed below."))

					, column(width=12, align="left", style="font-family: Merriweather;"
				 		, h4("Cluster 1"))
					, column(width=12, align="center"
						, plotOutput(outputId="cluster_1_plot", width=1800, height=1350))
					, column(width=12, align="left", style="font-family: Raleway; padding: 20px"
						, DT::dataTableOutput(outputId="cluster_1_dt", width="auto", height="auto"))


					, column(width=12, align="left", style="font-family: Merriweather;"
				 		, h4("Cluster 2"))
					, column(width=12, align="center"
						, plotOutput(outputId="cluster_2_plot", width=1800, height=1350))
					, column(width=12, align="left", style="font-family: Raleway; padding: 20px"
						, DT::dataTableOutput(outputId="cluster_2_dt", width="auto", height="auto"))


					, column(width=12, align="left", style="font-family: Merriweather;"
				 		, h4("Cluster 3"))
					, column(width=12, align="center"
						, plotOutput(outputId="cluster_3_plot", width=1800, height=1350))
					, column(width=12, align="left", style="font-family: Raleway; padding: 20px"
						, DT::dataTableOutput(outputId="cluster_3_dt", width="auto", height="auto"))


					, column(width=12, align="left", style="font-family: Merriweather;"
				 		, h4("Cluster 4"))
					, column(width=12, align="center"
						, plotOutput(outputId="cluster_4_plot", width=1800, height=1350))
					, column(width=12, align="left", style="font-family: Raleway; padding: 20px"
						, DT::dataTableOutput(outputId="cluster_4_dt", width="auto", height="auto"))


					, column(width=12, align="left", style="font-family: Merriweather;"
				 		, h4("Cluster 5"))
					, column(width=12, align="center"
						, plotOutput(outputId="cluster_5_plot", width=1800, height=1350))
					, column(width=12, align="left", style="font-family: Raleway; padding: 20px"
						, DT::dataTableOutput(outputId="cluster_5_dt", width="auto", height="auto"))
				)
			)



			, tabPanel ("Data"  # Contains the gene scores of every gene in the network, as well as an accompanying scatter plot
				, fluidRow (
					column(width=12, align="center", style="font-family: Merriweather; padding: 20px;"
						, h2("Gene Dataset Scores"))
					, column(width=12, align="center", style="font-family: Raleway;"
						, plotOutput(outputId="scores_plot", width=1280, height=720))
					, column(width=12, align="center", style="font-family: Raleway; padding: 20px;"
						, DT::dataTableOutput(outputId="genes_scores_dt", width="auto", height="auto"))
				)
			)

			, tabPanel ("Export"  # Allows for the export of data: graphics and data frames
				, fluidRow (
					column(width=12, align="left", style="font-family: Merriweather; padding: 20px;"
						, h3("Please select the file(s) you would like to export below:"))
					, column(width=2, align="left", style="font-family: Raleway;"
						, selectInput(inputId="export.select", label=NULL, choices=c("Gene List", "Gene Scores", "Gene Interactions")))
					, column(width=10, align="left", style="font-family: Raleway;"
						, downloadButton(outputId="export.button", label="Export"))
				)
			)

		)
	)
)

server <- function(input, output) {
	# Title panel will be hidden after either the Enter or Import buttons are pressed
	output$hide_title <- eventReactive(eventExpr=input$search_term.button, valueExpr=TRUE, ignoreInit = TRUE)
	output$hide_title <- eventReactive(eventExpr=input$import.button, valueExpr=TRUE, ignoreInit = TRUE)
	outputOptions(x=output, name="hide_title", suspendWhenHidden=FALSE)
	
	observeEvent (input$search_term.button | input$import.button, {
		if (input$search_term.button==0 && input$import.button==0) { return() }  # This ensures that the event is not observed if no buttons are pressed - It is a necessary check and the program will crash if not present

		toggle(id="title_animation", anim=TRUE, animType="fade", time=0.5, condition=NULL)  # Fade animation
		showModal(modalDialog(p(""), title = "Your request is being processed...", footer=NULL, size="l", easyClose=FALSE))  # Pop-up dialog to prevent user from clicking buttons again

		valid_input <- FALSE

		if (input$search_term.button == 1) {  # The enter button was pressed with a text query
			print("User entered a search term...")
			term <<- input$search_term.text
			initialize(term)
			if (term != "") {
				valid_input = TRUE
				genes <- read.csv(file="data/cross_checked_genes.csv", header=FALSE)
			}
		}
		else if (input$import.button == 1) {  # The import button was pressed and the user attached their own gene dataset
			print("User imported a set of genes...")
			in_file <- input$import.file
			if (!is.null(in_file)) {
				valid_input = TRUE
				genes <- read.csv(file=in_file$datapath, header=FALSE)  # Consider changing datapath to not be in a random temp folder
			}
			else {
				print("The imported file is empty of invalid.")
			}
		}

		if (valid_input) {
			# Initialization
			names(genes)[1] <- "gene"  # Rename data frame column to "gene"
			genes_mapped <- string_db$map(genes, "gene", removeUnmappedRows=TRUE)
			hits <- genes_mapped$STRING_id[1:nrow(genes_mapped)]
#			hits <- expand_network(hits)  # DO NOT USE YET! The expansion algorithm tends to add housekeeping genes to our network, which is not ideal.
			genes_interactions <<- string_db$get_interactions(hits)
			genes_scores <<- process_data(genes_interactions, hits)
			genes_interactions$from <<- map_names(genes_interactions$from)
			genes_interactions$to <<- map_names(genes_interactions$to)
			genes_interactions <<- genes_interactions[c(1, 2, 3, 5, 6, 8, 10, 12, 14, 16)]

			#---------- Leader Genes Tab ----------#
			leader_genes <- genes_scores[1:10,]
			output$leader_genes_dt <- DT::renderDataTable ({
				names(leader_genes)[names(leader_genes) == "gene"] <- "Gene Symbol"
				names(leader_genes)[names(leader_genes) == "STRING_id"] <- "STRING Identifier"
				names(leader_genes)[names(leader_genes) == "total_score"] <- "Global Score"
				DT::datatable(leader_genes) 
			})

			#---------- Network Tab ---------#
			output$network_plot <- renderPlot (
				{ string_db$plot_network(string_id=genes_scores$STRING_id, payload_id=NULL, required_score=950, add_link=FALSE, add_summary=FALSE) }
				, width=1800, height=1350)
			output$genes_interactions_dt <- DT::renderDataTable ({
				names(genes_interactions)[names(genes_interactions) == "from"] <- "Gene A"
				names(genes_interactions)[names(genes_interactions) == "to"] <- "Gene B"
				names(genes_interactions)[names(genes_interactions) == "neighborhood"] <- "Neighborhood"
				names(genes_interactions)[names(genes_interactions) == "fusion"] <- "Fusion"
				names(genes_interactions)[names(genes_interactions) == "cooccurence"] <- "Co-occurence"
				names(genes_interactions)[names(genes_interactions) == "coexpression"] <- "Co-expression"
				names(genes_interactions)[names(genes_interactions) == "experiments"] <- "Experimental"
				names(genes_interactions)[names(genes_interactions) == "database"] <- "Database"
				names(genes_interactions)[names(genes_interactions) == "textmining"] <- "Textmining"
				names(genes_interactions)[names(genes_interactions) == "combined_score"] <- "Combined Score"
				DT::datatable(genes_interactions, options=list(scrollX=T))
			})

			#---------- Clusters Tab ---------#
			clusters <- get_clusters(genes_scores$STRING_id)
			cluster_1 <- process_data(string_db$get_interactions(clusters[[1]]), clusters[1])
			cluster_2 <- process_data(string_db$get_interactions(clusters[[2]]), clusters[1])
			cluster_3 <- process_data(string_db$get_interactions(clusters[[3]]), clusters[1])
			cluster_4 <- process_data(string_db$get_interactions(clusters[[4]]), clusters[1])
			cluster_5 <- process_data(string_db$get_interactions(clusters[[5]]), clusters[1])

			output$cluster_1_plot <- renderPlot (
				{ string_db$plot_network(string_id=cluster_1$STRING_id, payload_id=NULL, required_score=950, add_link=FALSE, add_summary=FALSE) }
				, width=1800, height=1350)
			output$cluster_1_dt <- DT::renderDataTable (
				{ names(cluster_1)[names(cluster_1) == "gene"] <- "Gene Symbol"
				names(cluster_1)[names(cluster_1) == "STRING_id"] <- "STRING Identifier"
				names(cluster_1)[names(cluster_1) == "total_score"] <- "Global Score"
				DT::datatable(cluster_1) }
			)

			output$cluster_2_plot <- renderPlot (
				{ string_db$plot_network(string_id=cluster_2$STRING_id, payload_id=NULL, required_score=950, add_link=FALSE, add_summary=FALSE) }
				, width=1800, height=1350)
			output$cluster_2_dt <- DT::renderDataTable (
				{ names(cluster_2)[names(cluster_2) == "gene"] <- "Gene Symbol"
				names(cluster_2)[names(cluster_2) == "STRING_id"] <- "STRING Identifier"
				names(cluster_2)[names(cluster_2) == "total_score"] <- "Global Score"
				DT::datatable(cluster_2) }
			)

			output$cluster_3_plot <- renderPlot (
				{ string_db$plot_network(string_id=cluster_3$STRING_id, payload_id=NULL, required_score=950, add_link=FALSE, add_summary=FALSE) }
				, width=1800, height=1350)
			output$cluster_3_dt <- DT::renderDataTable (
				{ names(cluster_3)[names(cluster_3) == "gene"] <- "Gene Symbol"
				names(cluster_3)[names(cluster_3) == "STRING_id"] <- "STRING Identifier"
				names(cluster_3)[names(cluster_3) == "total_score"] <- "Global Score"
				DT::datatable(cluster_3) }
			)


			output$cluster_4_plot <- renderPlot (
				{ string_db$plot_network(string_id=cluster_4$STRING_id, payload_id=NULL, required_score=950, add_link=FALSE, add_summary=FALSE) }
				, width=1800, height=1350)
			output$cluster_4_dt <- DT::renderDataTable (
				{ names(cluster_4)[names(cluster_4) == "gene"] <- "Gene Symbol"
				names(cluster_4)[names(cluster_4) == "STRING_id"] <- "STRING Identifier"
				names(cluster_4)[names(cluster_4) == "total_score"] <- "Global Score"
				DT::datatable(cluster_4) }
			)


			output$cluster_5_plot <- renderPlot (
				{ string_db$plot_network(string_id=cluster_5$STRING_id, payload_id=NULL, required_score=950, add_link=FALSE, add_summary=FALSE) }
				, width=1800, height=1350)
			output$cluster_5_dt <- DT::renderDataTable (
				{ names(cluster_5)[names(cluster_5) == "gene"] <- "Gene Symbol"
				names(cluster_5)[names(cluster_5) == "STRING_id"] <- "STRING Identifier"
				names(cluster_5)[names(cluster_5) == "total_score"] <- "Global Score"
				DT::datatable(cluster_5) }
			)

			#---------- Data Tab ---------#
			output$scores_plot <- renderPlot (
				{ ggplot(genes_scores, aes(x=as.numeric(row.names(genes_scores)), y=genes_scores$total_score)) +
				geom_point(size=3) +
				geom_point(data=leader_genes, aes(x=as.numeric(row.names(leader_genes)), y=leader_genes$total_score), color="red", size=3) +
				xlab("Gene Index") +
				ylab("Global Score") +
				theme(axis.text=element_text(size=12), axis.title=element_text(size=16,face="bold")) }
				, width=1280, height=720
			)
			output$genes_scores_dt <- DT::renderDataTable (
				{ names(genes_scores)[names(genes_scores) == "gene"] <- "Gene Symbol"
				names(genes_scores)[names(genes_scores) == "STRING_id"] <- "STRING Identifier"
				names(genes_scores)[names(genes_scores) == "total_score"] <- "Global Score"
				DT::datatable(genes_scores) }
			)

			#---------- Export Tab ---------#
			data_to_export <- reactive ({
				switch (input$export.select
					, "Gene List"=genes_scores$gene
					, "Gene Scores"=genes_scores
					, "Gene Interactions"=genes_interactions)
			})
			output$export.button <- downloadHandler(
				filename=function() {
					if (input$export.select == "Gene List") { paste("genes_list", ".csv", sep="") }
					else if (input$export.select == "Gene Scores") { paste("genes_scores", ".csv", sep="") }
					else if (input$export.select == "Gene Interactions") { paste("genes_interactions", ".csv", sep="") }
				},
				content = function(file) {
					if (input$export.select == "Gene List") { write.table(data_to_export(), file, sep="\n", col.names=FALSE, row.names=FALSE) }
					else if (input$export.select == "Gene Scores") { write.csv(data_to_export(), file) }
					else if (input$export.select == "Gene Interactions") { write.csv(data_to_export(), file) }
				}
			)

			#  Pre-loads all navigation tabs so that users don't have to open the tab to begin loading data!
			#  However, this results in a prolonged white screen after processing the data, may need to look into this but it is a minor issue for now
			outputOptions(output, "network_plot", suspendWhenHidden = FALSE)
			outputOptions(output, "cluster_1_plot", suspendWhenHidden = FALSE)
			outputOptions(output, "cluster_2_plot", suspendWhenHidden = FALSE)
			outputOptions(output, "cluster_3_plot", suspendWhenHidden = FALSE)
			outputOptions(output, "cluster_4_plot", suspendWhenHidden = FALSE)
			outputOptions(output, "cluster_5_plot", suspendWhenHidden = FALSE)
			outputOptions(output, "scores_plot", suspendWhenHidden = FALSE)
		}
		else {
			print("The query is empty or invalid.")
		}

		removeModal()
	})
}

shinyApp(ui=ui, server=server)