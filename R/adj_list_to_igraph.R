#' adj_list_to_igraph
#' @description internal function that assembles all the KEGG data into a network/graph
#' @param adjl_G_S adjacency list containing relevant adjacent SNPs/KEGG genes
#' @return an igraph object, containing a network representing all the KEGG data
#'
adj_list_to_igraph <- function(adjl_G_S){
  ##creating IMON network by building up networks from SNP IMON level downwards
  ##creating edgelist dataframe for G_S (SNP -> gene):
  G_S <- utils::stack(adjl_G_S)
  G_S$values <- paste0("hsa: ", G_S$values) #adding "hsa:" identifier:
  G_S <- G_S[ , c(2,1)] #flipping columns
  ##creating edgelist dataframe for G_E (gene -> enzyme):
  G_E <- utils::stack(adjl_G_E)
  G_E <- G_E[, c(2,1)] #swap columns
  G_E$values <- gsub("\\[|\\]","",G_E$values) #ensure naming conventions match between levels
  ##find which genes that are mapped to snps are also mapped to enzymes
  G_E_clean <- G_E[which(G_E$ind %in% G_S$values),]
  ##creating edgelist dataframe for R_E (enzyme -> reaction):
  R_E <- utils::stack(adjl_R_E)
  R_E$values <- gsub("EC", "EC:", R_E$values) #ensuring naming conventions match between levels
  R_E$values <- gsub(" ", "", R_E$values) #ensuring naming conventions match between levels
  ##find which enzymes that are mapped to genes are also mapped to reactions
  R_E_clean <- R_E[which(R_E$values %in% G_E_clean$values),]
  ##creating edgelist dataframe for RP_R (reaction -> KEGG reaction pair):
  RP_R <- utils::stack(adjl_RP_R)
  ##find which reactions that are mapped to enzymes are also mapped to RPs
  RP_R_clean <- RP_R[which(RP_R$ind %in% R_E_clean$ind),]
  RP_R_clean <- RP_R_clean[ , c(2,1)] #flipping dataframe columns
  ##creating edgelist dataframe for RP_C (KEGG reaction pair <-> compound/metabolite):
  RP_C <- utils::stack(adjl_RP_C)
  ##find which RPs that are mapped to reactions are also mapped to compounds
  #RP_C_clean <- RP_C[which(RP_C$L1 %in% RP_R_clean$RP),]
  ##create igraph object for reaction pairs and compounds:
  g1 <- igraph::graph_from_data_frame(RP_C, directed = FALSE)
  g1 <- igraph::as.directed(g1, mode = "mutual") #add bidirection within the network
  ##create igraph object for reaction pairs and reactions
  g2 <- igraph::graph_from_data_frame(RP_R_clean, directed = FALSE)
  ##create igraph object for reactions and enzymes
  g3 <- igraph::graph_from_data_frame(R_E_clean, directed = TRUE)
  ##create igraph object for genes and enzymes
  g4 <- igraph::graph_from_data_frame(G_E_clean, directed = TRUE)
  ##create igraph object for genes and snps
  g5 <- igraph::graph_from_data_frame(G_S, directed = TRUE)
  ##merge all networks (RP_C and RP_R, R_E, E_G and G_S) into one single IMON:
  EL1 <- igraph::get.edgelist(g1) #get edgelists from networks
  EL2 <- igraph::get.edgelist(g2)
  EL3 <- igraph::get.edgelist(g3)
  EL4 <- igraph::get.edgelist(g4)
  EL5 <- igraph::get.edgelist(g5)
  ELU <- rbind(EL1, EL2, EL3, EL4, EL5) #take union of edgelists
  ELU1 <- ELU[!duplicated(ELU),] #take union of edgelists (again, remove any duplicates)
  ##creating overall IMON network as igraph object:
  G <- igraph::graph_from_edgelist(ELU1, directed =  TRUE)
  return(G)
}
