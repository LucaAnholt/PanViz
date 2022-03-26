#' adj_to_G
#'
#' @description Internal function that constructs an IMON (Integrated Multi-Omic Network) for an inputted adjacency list containing adjacency information between KEGG genes and queried SNPs.
#'
#' @param adjl_G_S - adjacency list containing relevant adjacencies between inputted SNPs and genes from KEGG
#'
#' @return igraph object representing total IMON for inputted SNPs
#'
#'
adjl_to_G <- function(adjl_G_S){
  G <- adj_list_to_igraph(adjl_G_S)
  G <- set_base_graph_attributes(G = G, colour_groups = FALSE)
  return(G)
}
