#' Create hashmap of associated KEGG IDs and compound/metabolite names
#'
#' @return hashmap/recursive list of associated KEGG IDs and compound/metabolite names
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
  pb <- tcltk::tkProgressBar(title = "Querying compound data from KEGG", min = 0, max = length(split_data), width = 400)
  progress <- function(n) tcltk::setTkProgressBar(pb, n, label=paste(round(n/length(split_data)*100,1),"% dowloaded"))
  opts <- list(progress = progress)
  Query_Compound_Data <- foreach::foreach(i = 1:length(split_data), .combine = 'c', .packages = c("tcltk"),.export = c("retry"),.options.snow = opts) %dopar% {
    library(KEGGREST)
    library(futile.logger)
    library(utils)
    PanViz:::retry(KEGGREST::keggGet(split_data[[i]]), maxErrors = 5, sleep = sleep)
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
  # saveRDS(names_list, "compound_names_hash.Rds")
  return(names_list)
}
