##dbSNP query function:
#' NCBI_dbSNP_query
#'
#' @param snp_list - list of SNPs to be queried via NCBI dbSNP API
#' @param progress_bar Boolean (default = TRUE) argument that controls whether or not a progress bar for calculations/KEGGREST API GET requests should be printed to the console
#'
#' @return - raw output from NCBI dbSNP API
NCBI_dbSNP_query <- function(snp_list, progress_bar){
  snp_list <- unique(gsub("rs", "", snp_list)) #remove rs from ID to query NCBI dbSNP
  split_data <- split(snp_list, ceiling(seq_along(snp_list)/100)) #split data for efficient query via NCBI Entrez API
  if(progress_bar == TRUE){
    pb <- utils::txtProgressBar(max = length(split_data), style = 3)
    raw_data_list <- list() #pre-initialise vector to capture raw dbSNP summary output
    for(i in seq_along(split_data)){ #query NCBI dbSNP API
      tryCatch(raw_data_list[[i]] <- rentrez::entrez_summary(db = "snp", id = split_data[[i]]), warning = function(w){})
      pctg <- paste(round(i/length(split_data) *100, 0), "% completed")
      utils::setTxtProgressBar(pb, i, label = pctg)
    }
    close(pb) #close progress bar
  }
  else{
    raw_data_list <- list() #pre-initialise vector to capture raw dbSNP summary output
    for(i in seq_along(split_data)){ #query NCBI dbSNP API
      tryCatch(raw_data_list[[i]] <- rentrez::entrez_summary(db = "snp", id = split_data[[i]]), warning = function(w){})
    }
  }
  ##removing recursion in list:
  raw_data <- unlist(raw_data_list, recursive = FALSE)
  return(raw_data)
}



