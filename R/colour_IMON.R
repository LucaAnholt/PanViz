#' colour network by categorical group levels
#'
#' @param G - igraph object containing uncoloured IMON
#'
#' @return - igraph object containing coloured IMON
#' @param progress_bar Boolean (default = TRUE) argument that controls whether or not a progress bar for calculations/KEGGREST API GET requests should be printed to the console
#'
colour_IMON <- function(G, progress_bar){
  snp_index <- igraph::V(G)[grepl("SNP", igraph::V(G)$type)] #get graph indexes for snps
  snp_colours <- igraph::V(G)[snp_index]$col #get colours for each snp (coloured by selected categorical variable)
  #get graph indexes for reactions:
  reaction_index <- which(igraph::V(G)$type %in% c("REACTION"))
  ##find index of first layer of RPs i.e. those directly connected to reactions:
  reaction_RP_edges <- igraph::E(G)[from(igraph::V(G)[reaction_index])]
  RPs <- unique(igraph::head_of(G, reaction_RP_edges)$name)
  RPs_index <- as.numeric(igraph::V(G)[RPs])
  directed_index <- which(igraph::V(G)$type %in% c("GENE", "ENZYME", "REACTION")) #get index of directed
  ##add indexes of first layer of RPs (directly connected to reactions):
  directed_index <- c(directed_index, RPs_index)
  ##get asymmetrical distance matrix of m snps to n (directly connected) nodes (m x n):
  directed_distances <- igraph::distances(G, v = igraph::V(G)[snp_index], to = igraph::V(G)[directed_index], mode = "out", weights = NA)
  metabolome_index <- which(igraph::V(G)$type %in% c("RP", "METABOLITE")) #get index of bidirectional nodes
  ##remove indexes of first layer of RPs (those directly connected to reactions):
  metabolome_index <- metabolome_index[!metabolome_index %in% RPs_index]
  ##define function for getting the graph indexes of directly connected nodes to a given snp:
  get_reachable_nodes <- function(snp, d){
    return(as.numeric(igraph::V(G)[names(d[snp,][which(!is.infinite(unname(d[snp, ])))])]))
  }
  ##pre-initialise an empty recursive list that will eventually contain all relative node colours:
  colour_list <- rep(list(c()), length(igraph::V(G)))
  ##add snp colours to colour_list:
  for(i in seq_along(snp_colours)){
    index <- as.numeric(snp_index[i])
    colour_list[[index]] <- snp_colours[i]
  }
  #set node colours based on reachable snps to genes, enzymes and reactions based on directed paths:
  for(i in seq_along(snp_index)){
    list_index <- get_reachable_nodes(i, directed_distances)
    colour <- snp_colours[i]
    for(j in seq_along(list_index)){
      index <- list_index[j]
      colour_list[[index]] <- c(colour_list[[index]], colour)
    }
  }
  #if a node has multiple connections to different coloured snps, merge these colours at this node:
  for(i in as.numeric(igraph::V(G)[directed_index])){
    if(length(colour_list[[i]]) > 1){
      colour_list[[i]] <- multi_hex_col_mix(colour_list[[i]])
    }
  }
  ##get asymmetrical distance matrix of m RPs (directly connected to reactions) to n metabolites/RPs (m x n):
  bidirected_distances <- igraph::distances(G, v = igraph::V(G)[RPs_index], to = igraph::V(G)[metabolome_index], mode = "all", weights = NA)
  ##function for finding distance from a given first layer RP to a metabolite
  get_reachable_node_distances <- function(snp, d){
    return(unname(d[snp, ][which(!is.infinite(unname(d[snp, ])))]))
  }
  ##function for standardising values between 0.1 and 0.9:
  ReScale <- function(x,first,last){(last-first)/(max(x)-min(x))*(x-min(x))+first}
  ##function for darkening hex colours:
  darken_colour <- function(cols, node_distances){
    #invisible(capture.output(cols1 <- readhex(file = textConnection(paste(cols, collapse = "\n")), class = "RGB"), file = "NUL"))
    cols1 <- colorspace::hex2RGB(cols)
    #transform to hue/lightness/saturation colorspace
    cols1 <- methods::as(cols1, "HLS")
    #multiplicative decrease of lightness
    cols1@coords[, "L"] <- as.vector(cols1@coords[, "L"]) * as.vector(node_distances)
    return(colorspace::hex(cols1))
  }
  ##get colours of first layer of RPs (directly connected to reactions/snps):
  RPs_colours <- unlist(colour_list[RPs_index])
  ##set metabolite node colours based on reachable RPs directlty connected to SNPs:
  for(i in seq_along(RPs_index)){
    list_index <- get_reachable_nodes(i, bidirected_distances) #get indexes
    node_distances <- get_reachable_node_distances(i, bidirected_distances) #get distances from RP (directly connected)
    node_distances <- ReScale(node_distances, 0.1, 0.9) #rescale distances between 0.01 and 0.3
    colour <- RPs_colours[i] #get colour of group
    list_of_colours <- rep(colour, length(node_distances)) #create vector of colour
    list_of_colours <- darken_colour(list_of_colours, node_distances) #edit colour by darkening by rescaled distance
    for(j in seq_along(list_index)){
      colour_edit <- list_of_colours[j]
      index <- list_index[j]
      colour_list[[index]] <- c(colour_list[[index]], colour_edit)
    }
  }
  if(progress_bar == TRUE){
    pb <- utils::txtProgressBar(max = length(metabolome_index), style = 3)
    ##merge all colours at every node in the metabolome:
    for(i in seq_along(as.numeric(igraph::V(G)[metabolome_index]))){
      index <- as.numeric(igraph::V(G)[metabolome_index])[i]
      if(length(colour_list[[index]]) > 1){
        colour_list[[index]] <- multi_hex_col_mix(colour_list[[index]])
      }
      pctg <- paste(round(i/length(metabolome_index) *100, 0), "% completed")
      utils::setTxtProgressBar(pb, i, label = pctg)
    }
    close(pb)
  }
  else if(progress_bar == FALSE){
    ##merge all colours at every node in the metabolome:
    for(i in seq_along(as.numeric(igraph::V(G)[metabolome_index]))){
      index <- as.numeric(igraph::V(G)[metabolome_index])[i]
      if(length(colour_list[[index]]) > 1){
        colour_list[[index]] <- multi_hex_col_mix(colour_list[[index]])
      }
    }
  }
  ##set node colour attributes:
  colour_vector <- unlist(colour_list)
  igraph::V(G)$color <- colour_vector
  if(progress_bar == TRUE){
    ##get out incident edge for nodes and colour by node colour:
    pb <- utils::txtProgressBar(max = length(igraph::V(G)), style = 3)
    for(i in seq_along(igraph::V(G))){
      node_colour <- igraph::V(G)[i]$col
      igraph::E(G)[from(igraph::V(G)[i])]$color <- node_colour
      pctg <- paste(round(i/length(igraph::V(G)) *100, 0), "% completed")
      utils::setTxtProgressBar(pb, i, label = pctg)
    }
    close(pb)
  }
  else if(progress_bar == FALSE){
    ##get out incident edge for nodes and colour by node colour:
    for(i in seq_along(igraph::V(G))){
      node_colour <- igraph::V(G)[i]$color
      igraph::E(G)[from(igraph::V(G)[i])]$color <- node_colour
    }
  }
  return(G) #return network with coloured attributes
}
