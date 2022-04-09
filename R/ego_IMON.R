#' ego_IMON
#' @description Internal function for trimming IMON to ego-centred (centred around SNPs) to specified order (pathway length from SNPs)
#' @param G - igraph object representing IMON
#' @param ego - the selected ego-centred path length
#' @return - ego-centred IMON set at desired path length
ego_IMON <- function(G, ego){
  ##get vertex indexes for inputted SNPs:
  snp_in = igraph::V(G)[which(igraph::V(G)$type %in% "SNP")]$name
  G <- igraph::make_ego_graph(graph = G, order = ego, nodes = igraph::V(G)[snp_in], mode = "out")
  ##merging list of egos into one single graph:
  subgraph_list_df <- lapply(G, igraph::as_data_frame) #create list of subgraph dataframes
  subgraph_df <- do.call(rbind, subgraph_list_df) #combine list of subgraph dataframes
  subgraph_df <- unique(subgraph_df) #ensure there are now duplicate edges
  G <- igraph::graph_from_data_frame(subgraph_df , directed = TRUE) #combine back into one single network
  #setting graph attributes:
  pal <- RColorBrewer::brewer.pal(n = 8, "Dark2") #getting colour palette for node colouring
  ##SNP attributes:
  igraph::V(G)[grepl("rs", igraph::V(G)$name)]$type <- "SNP"
  igraph::V(G)[grepl("rs", igraph::V(G)$name)]$color <- pal[1]
  #metabolite attributes:
  igraph::V(G)[grepl("C\\d{5}", igraph::V(G)$name)]$color <- pal[2]
  igraph::V(G)[grepl("C\\d{5}", igraph::V(G)$name)]$type <- "METABOLITE"
  #reaction attributes:
  igraph::V(G)[grepl("R\\d{5}", igraph::V(G)$name)]$color <- pal[3]
  igraph::V(G)[grepl("R\\d{5}", igraph::V(G)$name)]$type <- "REACTION"
  #enzyme attributes:
  igraph::V(G)[grepl("EC:", igraph::V(G)$name)]$color <- pal[4]
  igraph::V(G)[grepl("EC:", igraph::V(G)$name)]$type <- "ENZYME"
  #gene attributes
  igraph::V(G)[grepl("hsa:", igraph::V(G)$name)]$color <- pal[7]
  igraph::V(G)[grepl("hsa:", igraph::V(G)$name)]$type <- "GENE"
  #KEGG reaction pair attributes
  igraph::V(G)[grepl("C\\d{5}_C\\d{5}", igraph::V(G)$name)]$color <- pal[8]
  igraph::V(G)[grepl("C\\d{5}_C\\d{5}", igraph::V(G)$name)]$type <- "RP"
  ##set IDs of nodes:
  compound_names <- unname(compound_names_hash[igraph::V(G)[grepl("METABOLITE", igraph::V(G)$type)]$name])
  igraph::V(G)[grepl("METABOLITE", igraph::V(G)$type)]$ID <- compound_names
  igraph::V(G)[grepl("SNP", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("SNP", igraph::V(G)$type)]$name
  igraph::V(G)[grepl("GENE", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("GENE", igraph::V(G)$type)]$name
  igraph::V(G)[grepl("ENZYME", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("ENZYME", igraph::V(G)$type)]$name
  igraph::V(G)[grepl("REACTION", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("REACTION", igraph::V(G)$type)]$name
  igraph::V(G)[grepl("RP", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("RP", igraph::V(G)$type)]$name
  return(G)
}
