##algorithm to decompose non-connected IMONs
#' IMON decomposer
#'
#' @param G - igraph object containing non-fully connected IMON
#'
#' @return - list of igraph objects containing fully connected IMONs
decompose_IMON <- function(G){
  d = igraph::decompose(G) #decompose graph
  for(i in seq_along(d)){
    len = length(igraph::V(d[[i]])[grepl("SNP", igraph::V(d[[i]])$type)]$name) #remove any decomposed subgraph that doesn't contains snps
    if(len == 0){
      d[[i]] <- NA
    }
  }
  d = d[!is.na(d)]
  return(d)
}
