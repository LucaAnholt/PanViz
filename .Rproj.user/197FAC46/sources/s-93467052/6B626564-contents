##NCBI query check helper function
#' dbSNP_query_check
#'
#' @param query - raw query data from NCBI dbSNP API
#'
#' @return - vector containing either 0 (denoting successful query) or NA (unsuccessful query)
dbSNP_query_check <- function(query){
  if(is.null(query$error)){ #if there is no error in query return 0
    return(0)
  }
  else{ #if there is an error in query, return NA
    return(NA)
  }
}
