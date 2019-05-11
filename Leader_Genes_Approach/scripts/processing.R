#rm(list=ls())
#getwd()
#setwd("C:/Users/Jose/Documents/Ramapo/Bioinformatics/FinalProject/Leader_Genes_Approach/scripts")

## Jose Carrillo - Data Processing and Expansion - Advanced Bioinformatics - Spring 2019

# After mining for a preliminary set of cross-checked genes, we use R's STRINGdb package to perform analysis, but first we must process the data into a useable form.
# Gene processing enforces the condition of limiting our dataset to 400 genes because STRINGdb cannot accept more than that number when plotting network interactions.
# To process...
# 	1) We order the genes by their interaction score
# 	2) We only take, at most, the top 400 genes

library("reticulate", lib.loc="library")
library("STRINGdb", lib.loc="library")

string_db <- STRINGdb$new(version="10", species=9606, score_threshold=990, input_directory="")
aliases <- string_db$get_aliases()

expanded_genes <- c()  # Vector of genes that have been checked for association when expanding the network with STRING

"%!in%" <- Negate("%in%")  # Negation of %in% for easy usage :)

#---------- PRELIMINARY MINING ----------#

initialize <- function(term) {
	py$get_mesh_terms(term)
	py$export_mesh_terms()

	py$seed_preliminary_gene_uids()
	py$export_preliminary_gene_uids()

	for (i in 1:length(py$preliminary_gene_uids)) {
		if (cross_check_gene(py$preliminary_gene_uids[[i]])) {
			py$cross_checked_genes <- append(py$cross_checked_genes, py$preliminary_gene_uids[[i]])
		}
	}
	py$export_cross_checked_genes()	#CHANGECHANGE
	print("Finished initializing...")
}

#---------- DATA PROCESSING ----------#

process_data <- function(genes_interactions, hits) {  # Accepts a vector of STRING identifiers and returns a data frame with mapped names and total scores in descending order
	hits_with_interactions <- character(length(hits))

	count <- 0
	for (i in 1:nrow(genes_interactions)) {  # Loop removes orphan genes: genes without any interactions in the network
		if (genes_interactions[i, 1] %!in% hits_with_interactions) {  # Columns 1 and 2 contain STRING identifiers of the interaction
			hits_with_interactions[count] <- genes_interactions[i, 1]
			count <- count + 1
		}

		if (genes_interactions[i, 2] %!in% hits_with_interactions) {
			hits_with_interactions[count] <- genes_interactions[i, 2]
			count <- count + 1
		}
	}

	hits_with_interactions <- unique(hits_with_interactions[hits_with_interactions!=""])  # Remove empty elements from hits vector
	total_score <- numeric(length(hits_with_interactions))
	total_score[1:length(hits_with_interactions)] <- 0

	genes_scores <- data.frame(hits_with_interactions, total_score)
	colnames(genes_scores)[1] <- "STRING_id"  # Rename column

	genes_scores <- calculate_scores(genes_scores, genes_interactions)
	genes_scores$gene <- map_names(genes_scores$STRING_id)  # Adds a new column containing gene names/symbols
	genes_scores <- genes_scores[, c(3,1,2)]  # Re-ordeering columns to gene, STRING_id, total_score

	print("Finished processing...")
	return(genes_scores)
}

calculate_scores <- function(genes_scores, genes_interactions) {
	for (i in 1:nrow(genes_interactions)) {
		id_A <- genes_interactions[i, 1]  # STRING_id of the first gene in the interaction
		id_B <- genes_interactions[i, 2]  # STRING_id of the second gene in the interaction

		# We add the scores of the interaction to the respective row in genes_scores, column 16 contains the combined score value
		genes_scores[which(genes_scores$STRING_id == id_A),][2] <- genes_scores[which(genes_scores$STRING_id == id_A),][2] + genes_interactions[i, 16]
		genes_scores[which(genes_scores$STRING_id == id_B),][2] <- genes_scores[which(genes_scores$STRING_id == id_B),][2] + genes_interactions[i, 16]
	}

	genes_scores  # Data frame containing total interaction scores of each gene

	# A lot of manipulation...
	genes_scores <- genes_scores[order(-genes_scores$total_score),]  # Orders by descending total_score
	STRING_id <- genes_scores[,1]  # Vector for STRING identifiers
	total_score <- genes_scores[,2]  # Vector for total scores
	genes_scores <- as.data.frame(STRING_id)  # Create a new data frame (so that the rownames are ordered!)
	genes_scores$total_score <- total_score  # ...and some sort of bug forces me to declare the second column separately
	genes_scores <- genes_scores[1:400,]  # We only take the top 400 genes because STRINGdb cannot handle more than that
	genes_scores <- na.omit(genes_scores)  # And if we have less than 400 genes, we remove the null rows

	print("Calculated global scores...")
	return(genes_scores)  # Done!
}

get_clusters <- function(hits) {
	return(string_db$get_clusters(hits, algorithm="fastgreedy"))
	print("Resolved clusters...")
}

map_names <- function(STRING_id) {  # Accepts a vector of STRING identifiers and returns a vector of gene names/symbols corresponding to the identifier
	gene = character(length(STRING_id))
	for (i in 1:length(STRING_id)) {
		gene[i] <- aliases[aliases$STRING_id == STRING_id[i],][1,2]
	}

	print("Mapped gene symbols to STRING IDs...")
	return(gene)
}

#---------- EXPAND/CONVERGE ----------#

# This algorithm expands the current STRING network with predicted neighbors
# While it works perfectly (technically), it is not actually implemented (there is a code commented out in the Shiny app)
# This is because the expansion tends to include housekeeping/reference genes, and they co-occur with the search term in PubMed literature!!! (very bad)
# So, when we cross-check the housekeeping/reference gene and there is literature that has a co-occurrence, it believes that there is a functional association but we know that this is a false positive...
# Needless, to say, I cannot use it for now
expand_network <- function(hits) {
	print("Expanding network...")
	neighbors <- string_db$get_neighbors(hits)
	neighbors <- neighbors[neighbors %!in% hits]  # We do not take neighbors that are already in our data frame
	gene <- map_names(neighbors)

	#py$preliminary_gene_uids <- fromJSON(file="../data/preliminary_gene_uids.json")

	gene <- gene[gene %!in% py$preliminary_gene_uids]  # We do not consider preliminary genes since we have already crosschecked them all
	gene <- gene[gene %!in% expanded_genes]  # We do not consider neighbor genes that we have already crosschecked

	expanded_genes <<- append(expanded_genes, gene)  # The remaining neighbor genes are appended to expanded_genes, so that we do not crosscheck them again
	
	if(length(gene) == 0) {
		return(hits)
	}

	new_genes <- c()
	for (i in 1:length(gene)) {
		if (py$crosscheck_gene(gene[i]) == TRUE) {
			new_genes <- append(new_genes, gene[i])
		}
	}

	if (length(new_genes) != 0) {
		hits <- append(hits, string_db$mp(new_genes))
		expand_network(hits)
	}
	else {
		"Convergence reached..."
		return(hits)
	}
}

#---------- DEBUG ---------#