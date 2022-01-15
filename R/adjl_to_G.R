##adjaceny lists to graph function:
#' adj_to_G
#'
#' @param adjl_G_S
#'
#' @return igraph object representing total IMON for inputted SNPs
adjl_to_G <- function(adjl_G_S){
  cat("Generating IMON - this might take a few moments \n")
  ##creating IMON network by building up networks from SNP IMON level downwards
  ##creating edgelist dataframe for G_S (SNP -> gene):
  G_S <- stack(adjl_G_S)
  G_S$values <- paste0("hsa: ", G_S$values) #adding "hsa:" identifier:
  G_S <- G_S[ , c(2,1)] #flipping columns
  ##creating edgelist dataframe for G_E (gene -> enzyme):
  G_E <- stack(adjl_G_E)
  G_E <- G_E[, c(2,1)] #swap columns
  G_E$values = gsub("\\[|\\]","",G_E$values) #ensure naming conventions match between levels
  ##find which genes that are mapped to snps are also mapped to enzymes
  G_E_clean <- G_E[which(G_E$ind %in% G_S$values),]
  ##creating edgelist dataframe for R_E (enzyme -> reaction):
  R_E <- stack(adjl_R_E)
  R_E$values <- gsub("EC", "EC:", R_E$values) #ensuring naming conventions match between levels
  R_E$values <- gsub(" ", "", R_E$values) #ensuring naming conventions match between levels
  ##find which enzymes that are mapped to genes are also mapped to reactions
  R_E_clean <- R_E[which(R_E$values %in% G_E_clean$values),]
  ##creating edgelist dataframe for RP_R (reaction -> KEGG reaction pair):
  RP_R <- stack(adjl_RP_R)
  ##find which reactions that are mapped to enzymes are also mapped to RPs
  RP_R_clean <- RP_R[which(RP_R$ind %in% R_E_clean$ind),]
  RP_R_clean <- RP_R_clean[ , c(2,1)] #flipping dataframe columns
  ##creating edgelist dataframe for RP_C (KEGG reaction pair <-> compound/metabolite):
  #RP_C <- reshape2::melt(adjl_RP_C, value.name = "Compounds")
  RP_C <- stack(adjl_RP_C)
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
  ELU1 <-  ELU[!duplicated(ELU),] #take union of edgelists (again, remove any duplicates)
  ##creating overall IMON network as igraph object:
  G <- igraph::graph_from_edgelist(ELU1, directed =  TRUE)
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
  compound_names <- unname(compound_names_hash[igraph::V(G)[grepl("METABOLITE", igraph::V(G)$type)]$name])
  igraph::V(G)[grepl("METABOLITE", igraph::V(G)$type)]$ID <- compound_names
  igraph::V(G)[grepl("SNP", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("SNP", igraph::V(G)$type)]$name
  igraph::V(G)[grepl("GENE", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("GENE", igraph::V(G)$type)]$name
  igraph::V(G)[grepl("ENZYME", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("ENZYME", igraph::V(G)$type)]$name
  igraph::V(G)[grepl("REACTION", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("REACTION", igraph::V(G)$type)]$name
  igraph::V(G)[grepl("RP", igraph::V(G)$type)]$ID <- igraph::V(G)[grepl("RP", igraph::V(G)$type)]$name
  return(G)
}
