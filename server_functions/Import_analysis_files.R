# Qplots Version 3.0
# Last modified on 2/05/2019
# Script: Import Data (server functions)
# Author: Naima BELMOKHTAR (Email: naima1503@gmail.com)
#*******************************************************

#########################################
# import requested functions
source("tools_functions.R")

 #########################################
 # Import Otu table
 #########################################
 otu_file <- reactive({
   if (is.null(input$otumat_input)) {
     return(NULL)
   } else {
     ## check whether a .txt file is uploaded
     validate(
       need(tools::file_ext(input$otumat_input$name) %in% c('txt'), 
            "Wrong File Format Uploaded")
     )	
     ## Import file
     otumat <- read.table(input$otumat_input$datapath, header=TRUE, comment.char = "", stringsAsFactors = F, sep = "\t", row.names = 1, check.names = FALSE)
     # Clean table from empty lines
     otumat <- otumat[!apply(is.na(otumat) | otumat =="",1,all),]
     ##order row && column names
     otumat <- otumat[order(rownames(otumat)), order(colnames(otumat))]
     
    return(otumat)
     }
 })
 
 #########################################
 # Import Taxanomy table
 #########################################
 taxa_file <- reactive({
   if (is.null(input$taxamat_input)) {
     return(NULL)
   } else {
     ## check whether a .CSV or .txt file is uploaded
     validate(
       need(tools::file_ext(input$taxamat_input$name) %in% c('tsv', 'txt'), 
            "Wrong File Format Uploaded")
     )	
     ## Import file
     taxamat <- read.table(input$taxamat_input$datapath, header=TRUE, comment.char = "", stringsAsFactors = F, sep = "\t", check.names = FALSE)
     # Clean table from empty lines
     taxamat <- taxamat[!apply(is.na(taxamat) | taxamat =="",1,all),]
     ##order row names
     taxamat <- taxamat[order(rownames(taxamat)), ]
     taxatab <- split.taxon(taxamat)
     taxatab <- taxatab[order(rownames(taxatab)),]
     taxatab
   }
 })
 
 #########################################
 # Import Meta table
 #########################################
 meta_file <- reactive({
   if (is.null(input$metamat_input)) {
     return(NULL)
   } else {
     ## check whether a .CSV or .txt file is uploaded
     validate(
       need(tools::file_ext(input$metamat_input$name) %in% c('txt'), 
            "Wrong File Format Uploaded")
     )	
     ## Import file
     metamat <- read.table(input$metamat_input$datapath, header=TRUE, comment.char = "", stringsAsFactors = F, sep = "\t", row.names = 1, check.names = FALSE)
     # Clean table from empty lines
     metamat <- metamat[!apply(is.na(metamat) | metamat =="",1,all),]
     #order row names
     metamat <- metamat[order(rownames(metamat)),]
     metamat_filtered <- metamat %>%  
       select(where(~ n_distinct(.) > 1)) %>%
       select(where(~ n_distinct(.) < length(rownames(metamat))))
     metamat_filtered
   }
 })
 
 #########################################
 # Import Tree table
 #########################################
 tree_file <- reactive({
   if (is.null(input$tree_file_input)) {
     return(NULL)
   } else {
     ## check whether a .CSV or .txt file is uploaded
     validate(
       need(tools::file_ext(input$tree_file_input$name) %in% c('nwk','NWK', 'tree', 'tre'), 
            "Wrong File Format Uploaded")
     )	
     ## Import file
     tree <- ape::read.tree(input$tree_file_input$datapath)
     #tree <- phyloseq::read_tree_greengenes(fixUploadedFilesNames(input$tree_file_input)$datapath)
     tree
     }
 })
 
 #########################################
 # Build phyloseq data
 #########################################
  physeq <- reactive({
   if (is.null(input$otumat_input) && is.null(input$taxamat_input) && is.null(input$metamat_input)) {
    #if (is.null(input$read_input)) {
     return(NULL)
   } else {
     OTU = otu_table(as.matrix(otu_file()), taxa_are_rows = T)
     TAX = tax_table(as.matrix(taxa_file()))
     SAM = sample_data(meta_file())
     
     if(!is.null(tree_file())){
       TRE = phy_tree(tree_file())
       physeq = phyloseq(OTU, TAX, SAM, TRE)
     } else physeq = phyloseq(OTU, TAX, SAM)
     physeq
   }
 })
 
  
 #########################################
 # View Input Table
 #########################################
 observe({
   if (is.null(input$otumat_input) || is.null(input$taxamat_input) || is.null(input$metamat_input)) {
     shinyjs::disable("read_input")
   } else {
     shinyjs::enable("read_input")
   }
 })
 
 observeEvent(input$read_input,{
                output$otu_tab_view <- renderDataTable(otu_file())
                output$taxa_tab_view <- renderDataTable(taxa_file())
                output$meta_tab_view <- renderDataTable(meta_file())
                  
              })
 
