##algorithm for colouring node/edges where colours are mixed for intercepting different coloured incident edges
#' colour IMON nodes and edges by chosen grouping variable
#'
#' @param G - igraph object containing IMON
#'
#' @return - igraph object containing IMON with coloured nodes/edge attributes by grouping variable
colour_network_by_groups <- function(G){
  #check if IMON has all possible bioentities:
  if(length(unique(V(G)$type)) != 6){
    return(G)
  }
  cat("Setting node and edge colours at SNP/Gene level\n")
  ##get snp -> gene edges:
  snp_gene_edges <- E(G)[from(V(G)[grepl("rs", V(G)$name)])]
  ##colouring snp -> gene edges based on snp vertex colours
  for(i in seq_along(snp_gene_edges)){
    colour <- tail_of(G, snp_gene_edges[[i]])$col
    E(G)[from(V(G)[grepl("rs", V(G)$name)])][i]$col <- colour
  }
  ##colouring gene vertices based on incident edge colours
  genes <- V(G)[grepl("hsa:", V(G)$name)]$name
  for(i in seq_along(genes)){
    gene <- genes[[i]]
    incident_snp_cols <- tail_of(G, incident_edges(G, V(G)[grepl(gene, V(G)$name)], mode = "in")[[1]])$col
    if(length(incident_snp_cols)>1){
      merged_col <- multi_hex_col_mix(incident_snp_cols)
      V(G)[grepl(gene, V(G)$name)]$col <- merged_col
    }
    else{
      V(G)[grepl(gene, V(G)$name)]$col <- incident_snp_cols
    }
  }
  cat("Setting node and edge colours at Gene/Enzyme level\n")
  ##get gene -> enzyme edges:
  gene_enzyme_edges <- E(G)[from(V(G)[grepl("hsa:", V(G)$name)])]
  ##colouring gene -> enzyme edges based on gene vertex colours
  for(i in seq_along(gene_enzyme_edges)){
    colour <- tail_of(G, gene_enzyme_edges[[i]])$col
    E(G)[from(V(G)[grepl("hsa:", V(G)$name)])][i]$col <- colour
  }
  ##colouring enzyme vertices based on incident edge colours
  enzymes <- V(G)[grepl("EC:", V(G)$name)]$name
  for(i in seq_along(enzymes)){
    enzyme <- enzymes[[i]]
    incident_gene_cols <- tail_of(G, incident_edges(G, V(G)[grepl(enzyme, V(G)$name)], mode = "in")[[1]])$col
    if(length(incident_gene_cols)>1){
      merged_col <- multi_hex_col_mix(incident_gene_cols)
      V(G)[grepl(enzyme, V(G)$name)]$col <- merged_col
    }
    else{
      V(G)[grepl(enzyme, V(G)$name)]$col <- incident_gene_cols
    }
  }
  cat("Setting node and edge colours at Enzyme/Reaction level\n")
  ##get enzyme -> reaction edges:
  enzyme_reaction_edges <- E(G)[from(V(G)[grepl("EC:", V(G)$name)])]
  ##colouring enzyme -> reaction edges based on enzyme vertex colours
  for(i in seq_along(enzyme_reaction_edges)){
    colour <- tail_of(G, enzyme_reaction_edges[[i]])$col
    E(G)[from(V(G)[grepl("EC:", V(G)$name)])][i]$col <- colour
  }
  ##colouring reaction vertices based on incident edge colours
  reactions <- V(G)[grepl("R\\d{5}", V(G)$name)]$name
  for(i in seq_along(reactions)){
    reaction <- reactions[[i]]
    incident_enzyme_cols <- tail_of(G, incident_edges(G, V(G)[grepl(reaction, V(G)$name)], mode = "in")[[1]])$col
    if(length(incident_enzyme_cols)>1){
      merged_col <- multi_hex_col_mix(incident_enzyme_cols)
      V(G)[grepl(reaction, V(G)$name)]$col <- merged_col
    }
    else{
      V(G)[grepl(reaction, V(G)$name)]$col <- incident_enzyme_cols
    }
  }
  cat("Setting node and edge colours at Reaction/Reaction Pair level\n")
  ##get reaction -> RP edges:
  reaction_RP_edges <- E(G)[from(V(G)[grepl("R\\d{5}", V(G)$name)])]
  ##colouring reation -> RP edges based on reaction vertex colours
  for(i in seq_along(reaction_RP_edges)){
    colour <- tail_of(G, reaction_RP_edges[[i]])$col
    E(G)[from(V(G)[grepl("R\\d{5}", V(G)$name)])][i]$col <- colour
  }
  ##colouring RP vertices (only those directly connected to reactions) based on incident edge colours
  RPs <- head_of(G, reaction_RP_edges)$name
  for(i in seq_along(RPs)){
    RP <- RPs[[i]]
    incident_RP_cols <- tail_of(G, incident_edges(G, V(G)[grepl(RP, V(G)$name)], mode = "in")[[1]])$col
    incident_RP_cols <- incident_RP_cols[!is.na(incident_RP_cols)]
    if(length(incident_enzyme_cols)>1){
      merged_col <- multi_hex_col_mix(incident_RP_cols)
      V(G)[grepl(RP, V(G)$name)]$col <- merged_col
    }
    else{
      V(G)[grepl(RP, V(G)$name)]$col <- incident_RP_cols
    }
  }
  cat("Setting node and edge colours at Reaction Pair/Metabolite level\n")
  ##get all unique node colours
  all_node_colours <- unique(V(G)$col)[!is.na(unique(V(G)$col))]
  ##mix these colours to get one unique colour for the metabolome:
  mixed_colours <- multi_hex_col_mix(all_node_colours)
  ##find uncoloured RPs:
  uc_RPs <- V(G)[grepl("RP", V(G)$type)][which(is.na(V(G)[grepl("RP", V(G)$type)]$col))]
  V(G)[uc_RPs]$col <- mixed_colours #colour these with mixed colour
  E(G)[from(V(G)[uc_RPs])]$col <- mixed_colours
  ##colour all metabolites as mixed network colours:
  V(G)[grepl("METABOLITE", V(G)$type)]$col <- mixed_colours
  E(G)[from(V(G)[grepl("METABOLITE", V(G)$type)])]$col <- mixed_colours
  cat("Done setting node and edge colours\n")
  return(G)
}
