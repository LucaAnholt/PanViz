##enzyme -> reactions (adjl_R_E) adjacency list helper function:
#' Title
#'
#' @param queried_data - queried kegg data
#'
#' @return vector containing relevant adjacent enzymes and reactions
#'
adj_R_E <- function(queried_data){
  enzyme <- queried_data$ENZYME
  enzyme <- paste0("EC ", enzyme, sep = "") #add enzyme commission identifier
  return(enzyme)
}
