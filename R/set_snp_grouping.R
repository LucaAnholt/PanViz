##algorithm for setting SNP grouping names/colours by chosen grouping categorical variable
#' snp grouping by chosen categorical variable
#'
#' @param G - igraph object containing IMON
#' @param unique_group_names - vector containing unique grouping variable names
#' @param unique_group_cols - vector containing unique grouping colours for each variable
#' @param group_snps - snps split by each variable/group
#'
#' @return - igraph object containing IMON with labelled and coloured snps by grouping variable
set_snp_grouping <- function(G, unique_group_names, unique_group_cols, group_snps){
  for(i in seq_along(unique_group_names)){
    igraph::V(G)[which(igraph::V(G)$name %in% group_snps[[unique_group_names[[i]]]])]$group <- unique_group_names[[i]]
    igraph::V(G)[which(igraph::V(G)$name %in% group_snps[[unique_group_names[[i]]]])]$col <- unique_group_cols[[i]]
  }
  snp_index <- igraph::V(G)[grepl("SNP", igraph::V(G)$type)] #get graph indexes for snps
  non_snp_index <- igraph::V(G)[-snp_index]
  dm <- igraph::distances(G, v = igraph::V(G)[snp_index], to = igraph::V(G)[non_snp_index], mode = "out", weights = NA)
  return(G)
}
