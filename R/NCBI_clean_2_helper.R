##NCBI gene clean function 2
#' NCBI_clean_2
#'
#' @param queried_data - rentrez object queried from NCBI
#'
#' @return return chromosome location, start and end position of gene from NCBI query
#'
NCBI_clean_2 <- function(queried_data){
  if(length(queried_data) == 0){ #ensure that every item has genomic info
    chr_n <- NA
    chr_start <- NA
    chr_stop <- NA
  }
  else{
    chr_n <- queried_data$chrloc #separate genomic info
    chr_start <- queried_data$chrstart
    chr_stop <- queried_data$chrstop
  }
  return(list(chr_n, chr_start, chr_stop))
}
