#' Get_Kegg_Data
#' @description This function pulls the latest version of the KEGG database
#' @param CPU The number of cores to use when making KEGGREST API Get requests (default = 2). If CPU > 1, parallel requests will be made.
#' @param sleep The amount of sleep (seconds) between a potential caught API error and the next attempt (default = 5)
#' @return Rds files containing adjacency lists for all KEGG network data
#' @export
#' @examples \dontrun{
#' Get_Kegg_Data()
#' }
Get_Kegg_Data <- function(CPU = c(2,1), sleep = 5){
  ##Querying KEGG for metabolite, reaction and enzyme data:
  data <- retry(Reactions_Get_All(CPU = CPU, sleep = sleep), maxErrors = 3, sleep = sleep)
  adjl_R_E <- data[[1]]
  adjl_RP_C <- data[[2]]
  adjl_RP_R <- data[[3]]
  ##Querying KEGG for gene data + query gene locations from NCBI
  data <- retry(Genes_Get_All(CPU = CPU, sleep = sleep), maxErrors = 3, sleep = sleep)
  adjl_G_E <- data[[1]]
  Gene_Locations <- data[[2]]
  ##Querying compound names as hash:
  compound_names_hash <- retry(get_compound_hashmap(CPU = CPU, sleep = sleep), maxErrors = 3, sleep = sleep)
  #save to internal sys data:
  usethis::use_data(adjl_G_E, adjl_R_E, adjl_RP_C, adjl_RP_R, Gene_Locations, compound_names_hash, internal = TRUE, overwrite = TRUE)
  cat("Done! KEGG database successfully updated")
}
