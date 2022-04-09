#' get IMON with SNP and or all network vertices coloured by group variables
#' (either studies or phenotypes)
#' @description This function constructs an IMON (Integrated Multi-Omic Network)
#' with SNPs/or whole network coloured by selected categorical levels (either
#' studies or phenotypes)
#' @param dataframe A dataframe including 3 columns in the following order and
#' with the following names: snps, studies, traits (all character vectors)
#' @param groupby Choose whether to group SNP and or network colouring by either
#' studies or traits
#' @param ego This dictates what length order ego-centred network should be
#' constructed. If set to 5 (default and recommended), an IMON with the first
#' layer of the connected metabolome will be returned. If set above 5, the
#' corresponding extra layer of the metabolome will be returned. If set to 0
#' (not recommended) the fully connected metabolome will be returned. Note, this
#' cannot be set between 0 and 5.
#' @param save_file Boolean (default = FALSE) argument that indicates whether or
#' not the user wants to save the graph as an exported file in their current
#' working directory
#' @param export_type  This dictates the network data structure saved in your
#' working directory. By default this outputs an igraph object, however, you can
#' choose to export and save an edge list, graphml or GML file.
#' @param directory If set to "choose" this argument allows the user to interactively
#' select the directory of their choice in which they wish to save the constructed
#' IMON, else the file will be saved to the working directory "wd" by default
#' @param colour_groups Boolean (default = FALSE) chooses whether or not to colour
#' the whole network by grouping variables
#' @param progress_bar Boolean (default = TRUE) argument that controls whether or
#' not a progress bar for calculations/KEGGREST API GET requests should be printed
#' to the console
#' @return An igraph object containing the constructed IMON with coloured SNPs/and
#' or whole network by selected grouping variable
#' @export
#'
#' @examples
#' ##getting GWAS Catalog association tsv file and cleaning up using
#' ##GWAS_catalog_tsv_to_dataframe function:
#' path <- system.file("extdata",
#'   "gwas-association-downloaded_2021-09-13-EFO_1000649.tsv",
#'    package="PanViz")
#' df <- PanViz::GWAS_data_reader(file = path,
#'   snp_col = "SNPS",
#'   study_col = "STUDY",
#'   trait_col = "DISEASE/TRAIT")
#' ##creating uncoloured IMON:
#' G <- PanViz::get_grouped_IMON(dataframe = df,
#'   groupby = "studies",
#'   ego = 5,
#'   save_file = FALSE,
#'   colour_groups = FALSE)
#' ##creating IMON where vertices/edges are coloured by the variable study:
#' G <- PanViz::get_grouped_IMON(dataframe = df,
#'   groupby = "studies",
#'   ego = 5,
#'   save_file = FALSE,
#'   colour_groups = TRUE)
#'
get_grouped_IMON <- function(dataframe, groupby = c("studies", "traits"), ego = 5, save_file = c(FALSE, TRUE), export_type = c("igraph", "edge_list", "graphml", "gml"), directory = c("wd", "choose"), colour_groups = c(FALSE,TRUE), progress_bar = c(TRUE,FALSE)){
  ##test dimensions of inputted dataframe:
  if(dim(dataframe)[2] != 3){
    stop("Unexpected number of columns in dataframe. Dataframe must only include 3 columns: snps, studies and traits (in that order and naming convention)")
  }
  ##test if given snp vector is a vector of strings:
  if(sum(as.integer(as.vector(vapply(X = dataframe$snps, FUN = class, FUN.VALUE = as.character(1), USE.NAMES = FALSE)) %in% "character")) != length(dataframe$snps)){
    stop("Unexpected SNP vector. Inputted SNP vector must be supplied as character class")
  }
  ##test if given snp vector contains correct dbSNP accession number naming convention:
  if(sum(as.integer(stringr::str_detect(dataframe$snps, "rs"))) != length(dataframe$snps)){
    stop('Unexpected SNP column in inputted dataframe. Dataframe SNP column must use standard NCBI dbSNP accession naming convention (e.g. "rs185345278")')
  }
  if(missing(save_file)){
    save_file <- FALSE
  }
  if(missing(export_type)){
    export_type <- "igraph"
  }
  if(missing(progress_bar)){
    progress_bar <- TRUE
  }
  if(missing(colour_groups)){
    if(!missing(groupby)){
      colour_groups <- TRUE
    }
    else{
      colour_groups <- FALSE
    }
  }
  if(ego < 5 & ego != 0){ ##reset ego if below 5 and not set to special case of zero
    ego <- 5
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
  if(ego == 0){ ##if ego is set to zero
    warning("Cannot colour or group network when ego is set to zero.")
    return(get_IMON(snp_list = dataframe$snps, ego = 0, save_file = save_file, export_type = export_type, directory = directory))
  }
  ##get the categorical group by names (either studies or traits):
  unique_group_names <- unique(dplyr::pull(dataframe, column_index))
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
    snps <- data.frame(dataframe)[data.frame(dataframe)[,column_index] %in% value,]$snps
    return(snps)
  }
  ##create hash recursive list for finding associated snps for each group
  group_snps <- lapply(unique_group_names, filter_snps_by_group, dataframe, groupby)
  names(group_snps) <- unique_group_names
  ##get overall list of snps to query via NCBI dbSNP:
  snp_list <- dataframe$snps
  raw_data <- NCBI_dbSNP_query(snp_list, progress_bar)
  ##checking if all SNPs have been successfully queried:
  if(is.null(raw_data)){
    stop("None of the supplied SNPs could be queried via dbSNP")
  }
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
  ##finding matches between chromosome numbers (in SNPs and genes) and creating adjacency list:
  adjl_G_S <- snp_gene_map(Gene_Locations, snp_loc)
  if(length(adjl_G_S) == 0){ ##catch if no snps can be mapped to genes in database
    stop("None of the provided SNPs could be mapped to genes provided by the Kyoto Encyclopedia of Genes and Genomes")
  }
  ##creating IMON from all adjacency lists
  G <- adjl_to_G_grouped(adjl_G_S, unique_group_names, unique_group_cols, group_snps, colour_groups, ego = ego, progress_bar)
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
  }
  return(G) #returning IMON as igraph object for R manipulation/visualizations
}
