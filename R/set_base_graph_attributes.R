#' set_base_graph_attributes
#'
#' @param G igraph object containing KEGG network
#' @param colour_groups logical - whether or not user has indicated on colouring
#' the network by categorical variable i.e. study or trait/phenotype (only available
#' via PanViz::get_grouped_IMON())
#'
#' @return igraph object with node attributes set
#'
set_base_graph_attributes <- function(G, colour_groups){
  ##setting graph attributes:
  if(colour_groups == FALSE){
    pal <- RColorBrewer::brewer.pal(n = 8, "Dark2") #getting colour palette for node colouring
    igraph::V(G)[grepl("rs", igraph::V(G)$name)]$color <- pal[1]
    igraph::V(G)[grepl("C\\d{5}", igraph::V(G)$name)]$color <- pal[2]
    igraph::V(G)[grepl("R\\d{5}", igraph::V(G)$name)]$color <- pal[3]
    igraph::V(G)[grepl("EC:", igraph::V(G)$name)]$color <- pal[4]
    igraph::V(G)[grepl("hsa:", igraph::V(G)$name)]$color <- pal[7]
    igraph::V(G)[grepl("C\\d{5}_C\\d{5}", igraph::V(G)$name)]$color <- pal[8]
  }
  ##SNP attributes:
  igraph::V(G)[grepl("rs", igraph::V(G)$name)]$type <- "SNP"
  #metabolite attributes:
  igraph::V(G)[grepl("C\\d{5}", igraph::V(G)$name)]$type <- "METABOLITE"
  #reaction attributes:
  igraph::V(G)[grepl("R\\d{5}", igraph::V(G)$name)]$type <- "REACTION"
  #enzyme attributes:
  igraph::V(G)[grepl("EC:", igraph::V(G)$name)]$type <- "ENZYME"
  #gene attributes
  igraph::V(G)[grepl("hsa:", igraph::V(G)$name)]$type <- "GENE"
  #KEGG reaction pair attributes
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
  return(G)
}
