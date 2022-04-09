#' decompose_IMON
#'
#' @description This function returns a list of fully connected IMONs from a
#' single parent unconnected IMON.
#'
#' @param G - igraph object containing non-fully connected IMON
#'
#' @return - list of igraph objects, where each index contains a fully connected
#' IMON
#' @export
#' @examples
#' data("er_snp_vector")
#' G <- PanViz::get_IMON(snp_list = er_snp_vector, ego = 5, save_file = FALSE)
#' G_list <- decompose_IMON(G)
#'
decompose_IMON <- function(G){
  ##check if user has provided IMON (igraph object)
  if(missing(G)){
    stop("Please provide an unnconnected graph (IMON)")
  }
  ##attempt decompose using igraph:
  d <- igraph::decompose(G) #decompose graph
  ##if the network is already fully connected, show error:
  if(length(d) == 1){
    stop("The IMON you provided is already fully connected! Cannot decompose!")
  }
  ##if the network is not fully connected, separate the networks into a list of fully connected networks that fit IMON definition:
  for(i in seq_along(d)){
    len <- length(igraph::V(d[[i]])[grepl("SNP", igraph::V(d[[i]])$type)]$name) #remove any decomposed subgraph that doesn't contains snps
    if(len == 0){
      d[[i]] <- NA
    }
  }
  ##remove any potential NA values:
  d <- d[!is.na(d)]
  return(d)
}
