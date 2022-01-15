##NCBI query cleaning helper function
#' NCBI_clean
#'
#' @param queried_data - input queried NCBI gene data
#' @return remove genes with no genomic information from NCBI query
NCBI_clean <- function(queried_data){
  data <- queried_data$genomicinfo #get genomic info dataframe
  if(length(data)==0){
    return(NA) #if there's no genomic info, return NA
  }
  else{
    return(data) #else return genomic info
  }
}
