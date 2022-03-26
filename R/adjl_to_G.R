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
  cat("Setting IMON attributes \n")
  ##setting graph attributes:
  pal <- RColorBrewer::brewer.pal(n = 8, "Dark2") #getting colour palette for node colouring
  ##SNP attributes:
  igraph::V(G)[grepl("rs", igraph::V(G)$name)]$type <- "SNP"
  igraph::V(G)[grepl("rs", igraph::V(G)$name)]$col <- pal[1]
  #metabolite attributes:
  igraph::V(G)[grepl("C\\d{5}", igraph::V(G)$name)]$col <- pal[2]
  igraph::V(G)[grepl("C\\d{5}", igraph::V(G)$name)]$type <- "METABOLITE"
  #reaction attributes:
  igraph::V(G)[grepl("R\\d{5}", igraph::V(G)$name)]$col <- pal[3]
  igraph::V(G)[grepl("R\\d{5}", igraph::V(G)$name)]$type <- "REACTION"
  #enzyme attributes:
  igraph::V(G)[grepl("EC:", igraph::V(G)$name)]$col = pal[4]
  igraph::V(G)[grepl("EC:", igraph::V(G)$name)]$type = "ENZYME"
  #gene attributes
  igraph::V(G)[grepl("hsa:", igraph::V(G)$name)]$col = pal[7]
  igraph::V(G)[grepl("hsa:", igraph::V(G)$name)]$type = "GENE"
  #KEGG reaction pair attributes
  igraph::V(G)[grepl("C\\d{5}_C\\d{5}", igraph::V(G)$name)]$col = pal[8]
  igraph::V(G)[grepl("C\\d{5}_C\\d{5}", igraph::V(G)$name)]$type = "RP"
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
  return(G)
}
