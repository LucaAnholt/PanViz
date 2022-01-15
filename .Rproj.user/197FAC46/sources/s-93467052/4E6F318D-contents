##reaction (RP_C) adjacency list helper function:
#' adj_RP_C
#'
#' @param name - vector of all KEGG reaction names
#' @return vector containing relevant adjacent reaction pairs and compounds/metabolites
#'
adj_RP_C <- function(name){
  adjl_RP_C <- unlist(stringr::str_extract_all(name, "C\\d{5}")) #separating given reaction pair into associated metabolites/compounds
  return(adjl_RP_C)
}
