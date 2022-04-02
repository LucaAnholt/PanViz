#' get_Kegg_data
#' @description This function pulls the latest version of the KEGG database
#' @param CPU The number of cores to use when making KEGGREST API Get requests
#' (default = 2). If CPU > 1, parallel requests will be made. Note, a maximum
#' of two workers can be used.
#' @param sleep The amount of sleep (seconds) between a potential caught API error
#' and the next attempt (default = 5)
#' @param progress_bar Boolean (default = TRUE) argument that controls whether or
#' not a progress bar for calculations/KEGGREST API GET requests should be printed
#' to the console
#' @return Rds files containing adjacency lists for all KEGG network data
#' @importFrom usethis use_data
#' @export
#' @examples
#' PanViz::get_Kegg_data(CPU = 1, sleep = 5, progress_bar = FALSE)
get_Kegg_data <- function(CPU = 1, sleep = 5, progress_bar = c(TRUE, FALSE)){
  if(CPU > 2){
    stop("Can only use a maximum of 2 workers")
  }
  ##Querying KEGG for metabolite, reaction and enzyme data:
  data <- retry(Reactions_Get_All(CPU = CPU, sleep = sleep, progress_bar), maxErrors = 3, sleep = sleep)
  adjl_R_E <- append(data[[1]], adjl_R_E)
  adjl_RP_C <- append(data[[2]], adjl_RP_C)
  adjl_RP_R <- append(data[[3]], adjl_RP_R)
  ##Querying compound names as hash:
  compound_names_hash <- get_compound_hashmap()
  ##save to internal sys data:
  usethis::use_data(adjl_G_E, adjl_R_E, adjl_RP_C, adjl_RP_R, Gene_Locations, compound_names_hash, internal = TRUE, overwrite = TRUE)
}
