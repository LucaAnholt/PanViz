#' #' ego_IMON
#' #' @description internal function for trimming IMON to ego-centred (centred around SNPs) to specified order (pathway length from SNPs)
#' #'
#' #' @param G - igraph object representing IMON
#' #'
#' #' @param ego - the selected ego-centred path length
#' #'
#' #' @return - ego-centred IMON set at desired path length
#' ego_IMON_2 <- function(G, ego, grouped){
#'   if(grouped == TRUE){
#'     ##get attributes of parent graph:
#'     vertex_attr_dataframe <- data.frame(name = igraph::vertex_attr(G, "name", index = igraph::V(G)),
#'                                         type = igraph::vertex_attr(G, "type", index = igraph::V(G)),
#'                                         col = igraph::vertex_attr(G, "col", index = igraph::V(G)),
#'                                         group = igraph::vertex_attr(G, "group", index = igraph::V(G)),
#'                                         ID = unlist(igraph::vertex_attr(G, "ID", index = igraph::V(G)))
#'     )
#'     ##get edge attributes:
#'     edge_attr_dataframe <- data.frame(ID = 1:length(igraph::E(G)),
#'                                       col = igraph::edge_attr(G, "col", igraph::E(G))
#'     )
#'   }
#'   else if(grouped == FALSE){
#'     ##get attributes of parent graph:
#'     vertex_attr_dataframe <- data.frame(name = igraph::vertex_attr(G, "name", index = igraph::V(G)),
#'                                         type = igraph::vertex_attr(G, "type", index = igraph::V(G)),
#'                                         ID = unlist(igraph::vertex_attr(G, "ID", index = igraph::V(G)))
#'     )
#'   }
#'   ##get vertex indexes for inputted SNPs:
#'   snp_in = igraph::V(G)[which(igraph::V(G)$type %in% "SNP")]$name
#'   G <- igraph::make_ego_graph(graph = G, order = ego, nodes = igraph::V(G)[snp_in], mode = "out")
#'   ##merging list of egos into one single graph:
#'   subgraph_list_df <- lapply(G, igraph::as_data_frame) #create list of subgraph dataframes
#'   subgraph_df <- do.call(rbind, subgraph_list_df) #combine list of subgraph dataframes
#'   subgraph_df <- unique(subgraph_df) #ensure there are now duplicate edges
#'   G <- igraph::graph_from_data_frame(subgraph_df , directed = TRUE) #combine back into one single network
#'   ##switch over attributes from parent graph to new ego subgraph
#'   node_index <- which(vertex_attr_dataframe$name %in% igraph::V(G)$name)
#'   #setting graph attributes:
#'   igraph::V(G)$type <- vertex_attr_dataframe[node_index, ]$type
#'   igraph::V(G)$ID <- vertex_attr_dataframe[node_index, ]$ID
#'   if(grouped == TRUE){
#'   edge_index <- which(edge_attr_dataframe$ID %in% igraph::E(G)$ID)
#'   igraph::V(G)$group <- vertex_attr_dataframe[node_index, ]$group
#'   igraph::V(G)$col <- vertex_attr_dataframe[node_index, ]$col
#'   igraph::E(G)$col <- edge_attr_dataframe[edge_index, ]$col
#'   }
#'   return(G)
#' }
