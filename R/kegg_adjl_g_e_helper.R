##genes -> enzymes (adjl_G_E) adjacency list helper function:
#' adj_G_E
#'
#' @param queried_data
#'
#' @return vector containing relevant adjacent genes and enzymes
#'
adj_G_E <- function(queried_data){
  adjl_G_E <- queried_data$enzyme
  return(adjl_G_E)
}
