#' dbSNP query clean up function
#'
#' @description Internal function clean up raw SNP data queried from NCBI dbSNP via Entrez API depending on whether or not it could be successfully queried
#'
#' @param query - raw dbSNP query object
#'
#' @return - dataframe of separate chromosome number, position and ID
#'
dbSNP_query_clean <- function(query){
  if(as.integer(is.null(query$chrpos)) == 0){ #checking if query is successful
    chrpos <- query$chrpos #get genomic location
    chr <- stringr::str_extract(chrpos, "^([A-Za-z]{1,2}|[0-9]{1,2})") #extract chromosome number
    pos <- stringr::str_extract(chrpos, ":\\d+") #get genomic location
    pos <- gsub(":", "",pos) #remove colon
    ID <- query$uid #get SNP ID
  }
  else{
    chr <- NA #if SNP could not be queried return NAs
    pos <- NA
    ID <- NA
  }
  result = data.frame(snp_id = ID, chr_n = chr, chr_pos = pos) #concatenate results as dataframe
  return(result) #return dataframe
}
