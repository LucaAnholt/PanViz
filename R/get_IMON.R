#' get_IMON
#' @description Internal function that constructs an IMON (Integrated Multi-Omic Network)
#' for an inputted vector of SNPs and exports an igraph file.
#' @param snp_list A vector of SNPs (strings/characters) using standard NCBI dbSNP
#' accession number naming convention (e.g. "rs185345278")
#' @param ego This dictates what length order ego-centred network should be constructed.
#' If set to 5 (default and recommended), an IMON with the first layer of the connected
#' metabolome will be returned. If set above 5, the corresponding extra layer of the
#' metabolome will be returned. If set to 0 (not recommended) the fully connected
#' metabolome will be returned. Note, this cannot be set between 0 and 5.
#' @param save_file Boolean (default = FALSE) argument that indicates whether or
#' not the user wants to save the graph as an exported file in their current working
#' directory
#' @param export_type This dictates the network data structure saved in the chosen
#' directory. By default this outputs an igraph object, however, you can choose to
#' export and save an edge list, graphml or GML file.
#' @param directory If set to "choose" this argument allows the user to interactively
#' select the directory of their choice in which they wish to save the constructed
#' IMON, else the file will be saved to the working directory "wd" by default
#' @param progress_bar Boolean (default = TRUE) argument that controls whether or
#' not a progress bar for calculations/KEGGREST API GET requests should be printed
#' to the console
#'
#' @return An igraph object containing the constructed IMON
#' @export
#'
#' @examples
#' ##getting vector of SNPs to query:
#' data("er_snp_vector")
#' ##build IMON using vector:
#' G <- PanViz::get_IMON(snp_list = er_snp_vector, ego = 5, save_file = FALSE)
#'
get_IMON <- function(snp_list, ego = 5, save_file = c(FALSE, TRUE), export_type = c("igraph", "edge_list", "graphml", "gml"), directory = c("wd", "choose"), progress_bar = c(TRUE, FALSE)){
  ##test if given snp vector is a vector of strings:
  if(sum(as.integer(as.vector(vapply(X = snp_list, FUN = class, FUN.VALUE = as.character(1), USE.NAMES = FALSE)) %in% "character")) != length(snp_list)){
    stop("Unexpected SNP vector. Inputted SNP vector must be supplied as character class")
  }
  ##test if given snp vector contains correct dbSNP accession number naming convention:
  if(sum(as.integer(stringr::str_detect(snp_list, "rs"))) != length(snp_list)){
    stop('Unexpected SNP vector. Inputted SNP vector must use standard NCBI dbSNP accession naming convention (e.g. "rs185345278")')
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
  if(missing(snp_list)){
    stop("Please input a vector characters representing SNPs using the standard NCBI dbSNP accession naming convention, e.g. c('rs185345278', 'rs101')")
  }
  if(as.integer(export_type %in% c("igraph", "edge_list", "graphml", "gml")) == 0){
    stop('Unexpected export data type. export_type must be one of the following: "igraph", "edge_list", "graphml", "gml"')
  }
  if(ego < 5 & ego != 0){ ##reset ego if below 5 and not set to special case of zero
    ego <- 5
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
  gene_loc <- Gene_Locations
  row.names(gene_loc) <- NULL #remove indexing in row names
  ##finding matches between chromosome numbers (in SNPs and genes) and creating adjacency list:
  adjl_G_S <- snp_gene_map(gene_loc, snp_loc)
  if(length(adjl_G_S) == 0){
    stop("None of the provided SNPs could be mapped to genes provided by the Kyoto Encyclopedia of Genes and Genomes")
  }
  ##creating IMON from all adjacency lists
  G <- adjl_to_G(adjl_G_S)
  ##if user selects ego-centred network:
  if(ego != 0){
    G <- ego_IMON(G, ego)
  }
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
