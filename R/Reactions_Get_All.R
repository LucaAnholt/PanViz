#' Reactions_Get_All
#' @description This function constructs adjacency lists for compounds, reactions and enzymes listed within the KEGG database
#' @param CPU The number of cores to use when making KEGGREST API Get requests (default = 2). If CPU > 1, parallel requests will be made.
#' @param sleep The amount of sleep between a potential caught API error and the next attempt (default = 5)
#' @importFrom foreach %dopar%
#' @return Rds files for all relevant adjacency lists
#'
Reactions_Get_All <- function(CPU = c(2,1), sleep = 5){
  time = proc.time() #time process
  cat("Querying metabolite, reaction and enzyme data from KEGG\n")
  ##pulling reaction IDs from KEGG:
  Reaction_Raw_IDs <- KEGGREST::keggList("reaction")
  ##cleaning up raw data:
  Clean_Reactions <- attr(Reaction_Raw_IDs, "names")
  ##using raw reaction IDs to query KEGG:
  Query_Reaction_Data <- c()
  ##splitting data into chunks of 10 (max KEGG API search)
  split_data <- split(Clean_Reactions, ceiling(seq_along(Clean_Reactions)/10))
  ##Making parallel requests to KEGGREST API:
  cluster = parallel::makeCluster(CPU) #creating clusters
  doSNOW::registerDoSNOW(cluster)
  pb <- utils::txtProgressBar(max = length(split_data), style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  i <- NULL
  Query_Reaction_Data <- foreach::foreach(i = seq_along(split_data), .combine = 'c',.export = c("retry"),.options.snow = opts) %dopar% {
    retry(KEGGREST::keggGet(split_data[[i]]), maxErrors = 5, sleep = sleep)
  }
  parallel::stopCluster(cluster)
  close(pb)
  ##Cleaning up queried data and separating reaction pair data from compound data::
  Query_Reaction_Data <- lapply(Query_Reaction_Data, reaction_cleanup)
  ##Unlisting reaction pair IDs:
  name_list <- unique(unlist(lapply(Query_Reaction_Data, get_entry <- function(queried_data){if(all(!is.na(queried_data$RP))){return(queried_data$RP)}})))
  ##Creating adjacency list with reaction pairs and their related compounds
  adjl_RP_C <- lapply(name_list, adj_RP_C)
  names(adjl_RP_C) <- name_list
  ##making adjacency list with reactions and their related enzymes:
  ##getting vector of reaction names:
  reaction_names <- unique(unlist(lapply(Query_Reaction_Data, get_entry <- function(queried_data){return(queried_data$ENTRY)})))
  ##creating adjacency list for reactions -> reaction pairs (RP):
  adjl_RP_R <- lapply(Query_Reaction_Data, adj_RP_R)
  ##applying reaction names to adjl:
  names(adjl_RP_R) <- reaction_names
  ##removing any reactions with no associated reaction pairs:
  adjl_RP_R <- adjl_RP_R[!vapply(X = adjl_RP_R, FUN = function(x) all(is.na(x)), FUN.VALUE = as.logical(1))]
  ##creating adjacency list for enzymes -> reactions:
  adjl_R_E <- lapply(Query_Reaction_Data, adj_R_E)
  ##applying reaction names to adjl:
  names(adjl_R_E) <- reaction_names
  ##removing any reactions with no associated enzymes:
  adjl_R_E <- adjl_R_E[!vapply(X = adjl_RP_E, FUN = function(x) all(is.na(x)), FUN.VALUE = as.logical(1))]
  ##saving adjacency lists to selected directory:
  cat("Metabolite, reaction and enzyme adjacencies successfully queried - time elsapsed: ", (proc.time() - time)[[3]]/60, " minutes")
  cat("\n")
  return(list(adjl_R_E, adjl_RP_C, adjl_RP_R))
}



