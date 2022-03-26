#' Genes_Get_All
#' @description This function constructs adjacency lists for genes and enzymes listed within the KEGG database
#' @param CPU The number of cores to use when making KEGGREST API Get requests (default = 2). If CPU > 1, parallel requests will be made.
#' @param sleep The amount of sleep between a potential caught API error and the next attempt (default = 5).
#' @return Rds files for all relevant adjacency lists
#' @importFrom foreach %dopar%
Genes_Get_All <- function(CPU = c(2,1), sleep = 5){
  ##querying all human gene IDs from KEGG:
  all_Genes <- KEGGREST::keggList("hsa")
  ##Cleaning up the queried data:
  Clean_Genes <- attr(all_Genes, "names")
  ##splitting data into chunks of 10 (max KEGG API search)
  split_data <- split(Clean_Genes, ceiling(seq_along(Clean_Genes)/10))
  ##Querying all data from KEGG using the gene IDs:
  ##Making parallel requests to KEGGREST API:
  cluster = parallel::makeCluster(CPU) #creating clusters
  doSNOW::registerDoSNOW(cluster)
  pb <- utils::txtProgressBar(max = length(split_data), style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  i <- NULL
  Query_Genes_Data <- foreach::foreach(i = seq_along(split_data), .combine = 'c', .export = c("retry"),.options.snow = opts) %dopar% {
    Output <- c()
    retry(KEGGREST::keggGet(split_data[[i]]), maxErrors = 5, sleep = sleep)
  }
  parallel::stopCluster(cluster)
  close(pb)
  ##cleaning up this queried gene data from KEGG + filtering genes that have adjacent enzymes:
  Query_Genes_Data <- lapply(Query_Genes_Data, gene_cleanup)
  Query_Genes_Data <-  Query_Genes_Data[!is.na(Query_Genes_Data)]
  ##creating adjacency list for genes and enzymes:
  adjl_G_E <- lapply(Query_Genes_Data, adj_G_E)
  ##getting gene names:
  gene_ID <- unlist(lapply(Query_Genes_Data, get_name <- function(queried_data){return(paste0("hsa: ", queried_data$ENTRY[[1]]))}))
  ##adding gene names to adjacency list:
  names(adjl_G_E) <- gene_ID
  ##Finding the genomic locations of these genes:
  ##removing "hsa:" string:
  gene_ID <- gsub("hsa: ", "", gene_ID)
  ##Getting gene locations:
  Gene_Locations <- NCBI_Gene_Locations(Gene_Enterez_IDs = gene_ID)
  return(list(adjl_G_E, Gene_Locations))
}
