#' get IMON with SNP and or all network vertices coloured by group variables (either studies or phenotypes)
#' @description This function constructs an IMON (Integrated Multi-Omic Network) with SNPs/or whole network coloured by selected categorical levels (either studies or phenotypes)
#' @param dataframe A dataframe including 3 columns in the following order and with the following names: snps, studies, traits (all character vectors)
#' @param groupby Choose whether to group SNP and or network colouring by either studies or traits
#' @param ego This dictates what length order ego-centred network should be constructed. If set to 5 (default and recommended), an IMON with the first layer of the connected metabolome will be returned. If set above 5, the corresponding extra layer of the metabolome will be returned. If set to 0 (not recommended) the fully connected metabolome will be returned.
#' @param save_file Boolean (default = FALSE) argument that indicates whether or not the user wants to save the graph as an exported file in their current working directory
#' @param export_type  This dictates the network data structure saved in your working directory. By default this outputs an igraph object, however, you can choose to export and save an edge list, graphml or GML file.
#' @param directory If set to "choose" this argument allows the user to interactively select the directory of their choice in which they wish to save the constructed IMON, else the file will be saved to the working directory "wd" by default
#' @param colour_groups Boolean (default = FALSE) chooses whether or not to colour the whole network by grouping variables
#' @param progress_bar Boolean (default = TRUE) argument that controls whether or not a progress bar for calculations/KEGGREST API GET requests should be printed to the console
#' @return An igraph object containing the constructed IMON with coloured SNPs/and or whole network by selected grouping variable
#' @export
#'
get_grouped_IMON <- function(dataframe, groupby = c("studies", "traits"), ego = 5, save_file = c(FALSE, TRUE), export_type = c("igraph", "edge_list", "graphml", "gml"), directory = c("wd", "choose"), colour_groups = c(FALSE,TRUE), progress_bar = c(TRUE,FALSE)){
  ##test dimensions of inputted dataframe:
  if(dim(dataframe)[2] != 3){
    stop("Unexpected number of columns in dataframe. Dataframe must only include 3 columns: snps, studies and traits (in that order and naming convention)")
  }
  ##test if given snp vector is a vector of strings:
  if(sum(as.integer(as.vector(sapply(dataframe$snps, class)) %in% "character")) != length(dataframe$snps)){
    stop("Unexpected SNP vector. Inputted SNP vector must be supplied as character class")
  }
  ##test if given snp vector contains correct dbSNP accession number naming convention:
  if(sum(as.integer(stringr::str_detect(dataframe$snps, "rs"))) != length(dataframe$snps)){
    stop('Unexpected SNP column in inputted dataframe. Dataframe SNP column must use standard NCBI dbSNP accession naming convention (e.g. "rs185345278")')
  }
  if(missing(save_file)){
    save_file <- FALSE
  }
  ##test if given export_type argument is correct:
  if(save_file == FALSE){#ignore export_type argument if save_file is set to false
    invisible()
  }
  if(missing(export_type)){
    export_type <- "igraph"
  }
  if(missing(progress_bar)){
    progress_bar <- TRUE
  }
  if(as.integer(groupby %in% c("studies", "trais")) == 0){
    stop("Unexpected groupby value. This should be either 'studies' or 'traits'")
  }
  else if(as.integer(export_type %in% c("igraph", "edge_list", "graphml", "gml")) == 0){
    stop('Unexpected export data type. export_type must be one of the following: "igraph", "edge_list", "graphml", "gml"')
  }
  ##allow user to select which directory to save exported data file:
  if(missing(directory)){
    directory <- getwd()
  }
  if(directory == "wd"){
    directory <-  getwd()
  }
  if((directory == "choose") & (save_file == TRUE)){
    directory <- easycsv::choose_dir()
  }
  ##allow user to provide a custom file name:
  if(save_file == TRUE){
    filename <- readline(prompt = "Enter a custom file name, or press s to skip: ")
    if(filename == "s"){
      filename <-  "Exported_IMON"
    }
  }
  #get column index for selected categorical group:
  if(groupby == "studies"){
    column_index <- 2
  }
  if(groupby == "traits"){
    column_index <- 3
  }
  ##get the categorical group by names (either studies or traits):
  unique_group_names <- unique(as.character(dataframe[,column_index]))
  n <- length(unique_group_names) #get number of categories
  if(n > 50){
    stop('Exceeded the number of possible groups (50) - use PanViz::get_IMON() instead')
  }
  ##colour vector of 50 possible colours for grouping variables (colour blind friendly):
  col_vector <- c("#4b58a2","#9dbf46","#6640a3","#46ca79","#9c2e87","#4d9f3b","#d86cca",
              "#88c86b","#462875","#c1a830","#5982f0","#d6902c","#6f6cd3","#4d7b1f",
              "#b873d8","#4ec287","#d24283","#4ad3b6","#d43f61","#1fe1fb","#bd3439",
              "#2a9368","#d864b2","#7cbf76","#812862","#a2b356","#578ae2","#be6423",
              "#619ade","#db9c4c","#9b87da","#1f5113","#d694dc","#3e7e3c","#93579c",
              "#cdbc64","#952e56","#947d21","#db75a2","#747328","#ec605d","#8f541e",
              "#d4697b","#d4955a","#912641","#ce5a3b","#8c242e","#ce7058","#812614",
              "#d56467")
  unique_group_cols <- sample(col_vector, n) #create vector with associated colours for each categorical group
  ##function for filtering dataframe for snps associated with given categorical group
  filter_snps_by_group <- function(value, dataframe, groupby){
    snps <- dataframe[dataframe[,column_index] %in% value,]$snps
    return(snps)
  }
  ##create hash recursive list for finding associated snps for each group
  group_snps <- lapply(unique_group_names, filter_snps_by_group, dataframe, groupby)
  names(group_snps) <- unique_group_names
  ##get overall list of snps to query via NCBI dbSNP:
  snp_list <- dataframe$snps
  raw_data <- NCBI_dbSNP_query(snp_list, progress_bar)
  ##checking if all SNPs have been successfully queried:
  if(length(snp_list) > 1){ ##deal with vectorised input
    errors <- unname(unlist(lapply(raw_data, dbSNP_query_check)))
  }
  if(length(snp_list) == 1){ ##deal with non-vectorised input
    errors <- dbSNP_query_check(query = raw_data)
  }
  if(all(is.na(errors))){
    stop("None of the supplied SNPs could be queried via dbSNP")
  }
  ##clean up queried dbSNP data and return as dataframe with chromosome number, location and SNP ID
  snp_loc <- lapply(raw_data, dbSNP_query_clean)
  snp_loc <- do.call("rbind", snp_loc) #concatenate list of dataframes
  snp_loc <- stats::na.omit(snp_loc) #remove any rows with NAs
  snp_loc$chr_pos <- as.numeric(as.character(snp_loc$chr_pos))
  ##retrieving KEGG gene locations from sys data:
  row.names(Gene_Locations) <- NULL #remove indexing in row names
  cat("Mapping inputted SNPs to KEGG genes\n")
  ##finding matches between chromosome numbers (in SNPs and genes) and creating adjacency list:
  adjl_G_S <- snp_gene_map(Gene_Locations, snp_loc)
  if(length(adjl_G_S) == 0){ ##catch if no snps can be mapped to genes in database
    stop("None of the provided SNPs could be mapped to genes provided by the Kyoto Encyclopedia of Genes and Genomes")
  }
  ##creating IMON from all adjacency lists
  G <- adjl_to_G_grouped(adjl_G_S, unique_group_names, unique_group_cols, group_snps, colour_groups, ego = ego)
  ##saving graphs to desired output
  if(save_file == TRUE){
    if(export_type == "igraph"){
      saveRDS(G, file = paste0(directory, "/", filename, ".Rds"))
    }
    if(export_type == "edge_list"){
      igraph::write_graph(G, file = paste0(directory, "/", filename, ".txt"), format = "edgelist")
    }
    if(export_type == "graphml"){
      igraph::write_graph(G, file = paste0(directory, "/", filename, ".graphml"), format = "graphml")
    }
    if(export_type == "gml"){
      igraph::write_graph(G, file = paste0(directory, "/", filename, ".gml"), format = "gml")
    }
    cat(paste0("Done! IMON saved in ", directory, " as ", "'",filename,".",export_type,"'","\n"))
  }
  else{
    cat("Done!")
  }
  return(G) #returning IMON as igraph object for R manipulation/visualizations
}
