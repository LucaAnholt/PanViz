#' adjl_to_G_grouped
#' @description Internal function that constructs either a variable-coloured or uncoloured IMON (Integrated Multi-Omic Network) for an inputted adjacency list containing adjacency information between KEGG genes and queried SNPs.
#' @param adjl_G_S - adjacency list containing relevant adjacencies between inputted SNPs and genes from KEGG
#' @param unique_group_names - a list of the unique group/variable names in the provided GWAS Catalog association file
#' @param unique_group_cols - a list of unique colours for each unique group/variable in the provided GWAS Catalog association file
#' @param group_snps - a recursive list containing the lists of SNPs belonging to each unique group/variable in the provided GWAS Catalog association file
#' @param colour_groups - boolean: whether or not user has chosen to colour the network by the unique group/variables in the provided GWAS Catalog association file
#' @param ego - the egocentric order (centred around the SNPs in the network) in which to build the network i.e. pathlength from SNPs downwards towards the metabolome
#' @param progress_bar - boolean: whether or not user has decided to have a progress bar print to the console
#'
#' @return - an igraph object containing the IMON
#'
#'
adjl_to_G_grouped <- function(adjl_G_S, unique_group_names, unique_group_cols, group_snps, colour_groups, ego, progress_bar){
  G <- adj_list_to_igraph(adjl_G_S)
  ##setting graph attributes:
  G <- set_snp_grouping(G, unique_group_names, unique_group_cols, group_snps)
  G <- set_base_graph_attributes(G = G, colour_groups = colour_groups)
  ##create ego-centred IMON to required path length selected by user
  G <- ego_IMON(G, ego)
  ##check if network length is reasonable for colouring:
  if(length(igraph::V(G)) > 2000){
    stop("IMON is too large to colour whole network by categorical levels")
  }
  ##colour graph by user selection
  if(colour_groups == TRUE){
    G <- set_snp_grouping(G, unique_group_names, unique_group_cols, group_snps) #if is fully connected colour IMON
    tryCatch(G <- colour_IMON(G, progress_bar), warning = function(w){})
  }
  ##set edge ID attribute:
  igraph::E(G)$ID <- seq_along(igraph::E(G))
  return(G)
}
