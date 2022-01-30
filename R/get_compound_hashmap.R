#' get_compound_hashmap
#' @description This function creates a hashmap of associated KEGG IDs and compound/metabolite names
#' @param CPU The number of cores to use when making KEGGREST API Get requests (default = 2). If CPU > 1, parallel requests will be made.
#' @param sleep The amount of sleep (seconds) between a potential caught API error and the next attempt (default = 5)
#' @return hashmap/recursive list of associated KEGG IDs and compound/metabolite names
#'
#' @examples \dontrun{
#' compound_hashmap <- get_compound_hashmap()
#' }
get_compound_hashmap <- function(CPU = c(2,1), sleep = 5){
  compounds_Raw_IDs <- KEGGREST::keggList("compound")
  ##cleaning up raw data:
  Clean_Compounds <- attr(compounds_Raw_IDs, "names")
  ##using raw reaction IDs to query KEGG:
  Query_Compound_Data <- c()
  ##splitting data into chunks of 10 (max KEGG API search)
  split_data <- split(Clean_Compounds, ceiling(seq_along(Clean_Compounds)/10))
  ##Making parallel requests to KEGGREST API:
  cluster = parallel::makeCluster(CPU) #creating clusters
  doSNOW::registerDoSNOW(cluster)
  pb <- utils::txtProgressBar(max = length(split_data), style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  Query_Compound_Data <- foreach::foreach(i = 1:length(split_data), .combine = 'c', .export = c("retry"),.options.snow = opts) %dopar% {
    retry(KEGGREST::keggGet(split_data[[i]]), maxErrors = 5, sleep = sleep)
  }
  parallel::stopCluster(cluster)
  close(pb)
  ##define function for cleaning KEGG data
  cleaner <- function(recursive_list){
    return(recursive_list$NAME[[1]])
  }
  ##preinitialise hash map (recursive list) with associated KEGG metabolite/compound id with associated name
  names_list <- c()
  for(i in seq_along(Query_Compound_Data)){
    names_list[[i]] <- cleaner(Query_Compound_Data[[i]])
  }
  ##get cleaned KEGG metabolite/compound IDs:
  compounds_names <- unlist(stringr::str_extract_all(Clean_Compounds, "C\\d{5}"))
  ##name hash map/recursive list
  names(names_list) <- compounds_names
  ##save as RDS
  return(names_list)
}
