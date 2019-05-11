#------------------- Leader Genes Approach - Advanced Bioinformatics - Spring 2019 - Jose Carrillo -------------------#

This application aims to identify functional genes of a particular disease, phenotype, or biological process that has been well-documented in NCBI. Out of the set of genes, those predicted to have the most central or critical role in the network are termed "leader" genes.

The application is written in Python and R, and the front-end is implemented as a Shiny R application.


#-----INITIALIZATION-----#
Please include the following libraries:
--Python--
1. Bio.Entrez
2. csv
3. json
4. time.sleep
5. urllib3.exceptions.HTTPError
6. xml.etree.ElementTree

--R--
1. DT
2. ggplot2
3. reticulate
4. shiny
5. shinyjs
6. STRINGdb (from Bioconductor, I already included STRINGdb in the "library" subdirectory, so there is no need to install this!)

Please source the "app.R" script and run the following command:
	shinyApp(ui=ui, server=server)

Or, simply set your working directory to the location of "app.R" and run the script.


#-----INPUTTING DATA-----#
The title page gives you two options: 1) Enter your own search term and press "Enter" or 2) Import your own set of genes

Entering your own search term initializes the data mining script, which has a runtime of about 45 to 90 minutes! I have already ran some tests with the data mining script to generate gene lists that are ready to be imported. These can be found in the "data" subdirectory as .csv files. Importing these gene lists will decrease the runtime to about only 5 minutes or less, so please feel free to use them to test the front-end and analysis part of my application.

At any point that there is some data processing running, I have included printouts for nearly every processing function that tracks the progress. They should all be printed out on the R console.

Once the data has been mined (if you entered a search term) and processed, the screen will turn white for at most 2 minutes. This is due to the rendering of the STRING network plots, and it unfortunately also removes my "Your request is being processed..." prompt. Fear not, the application has neither crashed nor returned you a blank white screen, so please give it a bit of time!


#------VIEWING DATA-----#
The application contains five navigation tabs in total, four of which organize the data and the final tab is for exporting the data.

1. Leader Genes: Lists the top ten leader genes of the network as well as their global scores.

2. Network: Plots the STRING network and tabulates every single interaction in the network as well as the weight of each link.

3. Clusters: Plots the first five major clusters of the network and their respective interactions. As the genes are ordered by score, the leader genes of each cluster can also be determined.

4. Data: Plots the global scores of every gene, including the leader genes and tabulates this data. In the plot, the points in red mark the leader genes.

5. Export: Allows the user to export the list of genes used in the network, the table of all gene interactions, as well as the table of all gene scores.
