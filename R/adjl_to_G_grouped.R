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
  cat("Setting IMON attributes \n")
  ##setting graph attributes:
  ##SNP attributes:
  igraph::V(G)[grepl("rs", igraph::V(G)$name)]$type <- "SNP"
  ##metabolite attributes:
  igraph::V(G)[grepl("C\\d{5}", igraph::V(G)$name)]$type <- "METABOLITE"
  ##reaction attributes:
  igraph::V(G)[grepl("R\\d{5}", igraph::V(G)$name)]$type <- "REACTION"
  ##enzyme attributes:
  igraph::V(G)[grepl("EC:", igraph::V(G)$name)]$type <- "ENZYME"
  ##gene attributes
  igraph::V(G)[grepl("hsa:", igraph::V(G)$name)]$type <- "GENE"
  ##KEGG reaction pair attributes
  igraph::V(G)[grepl("C\\d{5}_C\\d{5}", igraph::V(G)$name)]$type <- "RP"
  ##set IDs of nodes:
  index <- which(names(compound_names_hash) %in% igraph::V(G)[grepl("METABOLITE", igraph::V(G)$type)]$name)
  index <- index[order(match(names(compound_names_hash[index]),igraph::V(G)[grepl("METABOLITE", igraph::V(G)$type)]$name))]
  compound_names <- unname(unlist(compound_names_hash[index]))
  igraph::V(G)[grepl("METABOLITE", igraph::V(G)$type)]$ID <- compound_names
  igraph::V(G)[grepl("SNP", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("SNP", igraph::V(G)$type)]$name
  igraph::V(G)[grepl("GENE", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("GENE", igraph::V(G)$type)]$name
  igraph::V(G)[grepl("ENZYME", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("ENZYME", igraph::V(G)$type)]$name
  igraph::V(G)[grepl("REACTION", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("REACTION", igraph::V(G)$type)]$name
  igraph::V(G)[grepl("RP", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("RP", igraph::V(G)$type)]$name
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
  else{
    G <- set_snp_grouping(G, unique_group_names, unique_group_cols, group_snps) #if is fully connected colour IMON
    pal <- RColorBrewer::brewer.pal(n = 8, "Dark2") #getting colour palette for node colouring
    igraph::V(G)[grepl("SNP", igraph::V(G)$type)]$col <- pal[1]
    igraph::V(G)[grepl("METABOLITE", igraph::V(G)$type)]$col <- pal[2]
    igraph::V(G)[grepl("REACTION", igraph::V(G)$type)]$col <- pal[3]
    igraph::V(G)[grepl("ENZYME", igraph::V(G)$type)]$col <- pal[4]
    igraph::V(G)[grepl("GENE", igraph::V(G)$type)]$col <- pal[7]
    igraph::V(G)[grepl("RP", igraph::V(G)$type)]$col <- pal[8]
    index <- which(names(compound_names_hash) %in% igraph::V(G)[grepl("METABOLITE", igraph::V(G)$type)]$name)
    index <- index[order(match(names(compound_names_hash[index]),igraph::V(G)[grepl("METABOLITE", igraph::V(G)$type)]$name))]
    compound_names <- unname(unlist(compound_names_hash[index]))
    igraph::V(G)[grepl("METABOLITE", igraph::V(G)$type)]$ID <- compound_names
    igraph::V(G)[grepl("SNP", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("SNP", igraph::V(G)$type)]$name
    igraph::V(G)[grepl("GENE", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("GENE", igraph::V(G)$type)]$name
    igraph::V(G)[grepl("ENZYME", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("ENZYME", igraph::V(G)$type)]$name
    igraph::V(G)[grepl("REACTION", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("REACTION", igraph::V(G)$type)]$name
    igraph::V(G)[grepl("RP", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("RP", igraph::V(G)$type)]$name
  }
  ##set edge ID attribute:
  igraph::E(G)$ID <- seq_along(igraph::E(G))
  return(G)
}
