##ego graph generating function 2
#' ego_IMON
#'
#' @param G - igraph object representing IMON
#' @param ego - the selected ego-centred path length
#'
#' @return - ego-centred IMON set at desired path length
ego_IMON_2 <- function(G, ego){
  ##get attributes of parent graph:
  vertex_attr_dataframe <- data.frame(name = vertex_attr(G, "name", index = V(G)),
                                         type = vertex_attr(G, "type", index = V(G)),
                                         col = vertex_attr(G, "col", index = V(G)),
                                         group = vertex_attr(G, "group", index = V(G)),
                                         ID = vertex_attr(G, "ID", index = V(G))
                                         )

  edge_attr_dataframe <- data.frame(ID = 1:length(E(G)),
                                      col = edge_attr(G, "col", E(G))
                                      )
  ##get vertex indexes for inputted SNPs:
  snp_in = V(G)[which(V(G)$type %in% "SNP")]$name
  G <- make_ego_graph(graph = G, order = ego, nodes = V(G)[snp_in], mode = "out")
  ##merging list of egos into one single graph:
  subgraph_list_df <- lapply(G, as_data_frame) #create list of subgraph dataframes
  subgraph_df <- do.call(rbind, subgraph_list_df) #combine list of subgraph dataframes
  subgraph_df <- unique(subgraph_df) #ensure there are now duplicate edges
  G <- graph_from_data_frame(subgraph_df , directed = TRUE) #combine back into one single network
  ##switch over attributes from parent graph to new ego subgraph
  node_index <- which(vertex_attr_dataframe$name %in% V(G)$name)
  edge_index <- which(edge_attr_dataframe$ID %in% E(G)$ID)
  #setting graph attributes:
  V(G)$type <- vertex_attr_dataframe[node_index, ]$type
  V(G)$col <- vertex_attr_dataframe[node_index, ]$col
  V(G)$group <- vertex_attr_dataframe[node_index, ]$group
  V(G)$ID <- vertex_attr_dataframe[node_index, ]$ID
  E(G)$col <- edge_attr_dataframe[edge_index, ]$col
  return(G)
}
