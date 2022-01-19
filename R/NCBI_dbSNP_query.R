##dbSNP query function:
#' NCBI_dbSNP_query
#'
#' @param snp_list - list of SNPs to be queried via NCBI dbSNP API
#'
#' @return - raw output from NCBI dbSNP API
NCBI_dbSNP_query <- function(snp_list){
  snp_list <- unique(gsub("rs", "", snp_list)) #remove rs from ID to query NCBI dbSNP
  split_data <- split(snp_list, ceiling(seq_along(snp_list)/100)) #split data for efficient query via NCBI Entrez API
  cat("Querying SNP locations from NCBI dbSNP \n")
  pb <- utils::txtProgressBar(max = length(split_data), style = 3)
  raw_data_list <- list() #pre-initialise vector to capture raw dbSNP summary output
  for(i in seq_along(split_data)){ #query NCBI dbSNP API ###api_key = "24f5a5e0cf387cfb54f726abc4458b62d109"
    raw_data_list[[i]] <- rentrez::entrez_summary(db = "snp", id = split_data[[i]])
    pctg <- paste(round(i/length(split_data) *100, 0), "% completed")
    utils::setTxtProgressBar(pb, i, label = pctg)
  }
  close(pb) #close progress bar
  ##removing recursion in list:
  raw_data <- unlist(raw_data_list, recursive = FALSE)
  return(raw_data)
}
