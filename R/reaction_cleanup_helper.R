##reaction KEGG data clean up helper function:
#' reaction_cleanup
#' @description This function helps to cleans up queried KEGG reaction recursive lists + separates compound/metabolite and reaction pair data into new sections
#' @param queried_data - input queried KEGG reaction data
#' @return Trimmed recursive lists containing queried KEGG reaction data
reaction_cleanup <- function(queried_data){
  ##deleting unnecessary data
  queried_data$DEFINITION <- NULL
  queried_data$PATHWAY <- NULL
  queried_data$MODULE <- NULL
  queried_data$ORTHOLOGY <- NULL
  queried_data$DBLINKS <- NULL
  ##reorganizing compound/RCLASS (reaction pair (RP)) data entries:
  if(length(queried_data$RCLASS) > 0){
    output <- gsub("RC\\d{5}\\s ", "", queried_data$RCLASS)
    output <- stringr::str_extract_all(output, "C\\d{5}")
    output <- as.vector(unlist(output))
    queried_data$COMPOUNDS <- output
    rp <- stringr::str_extract_all(queried_data$RCLASS, "C\\d{5}_C\\d{5}")
    queried_data$RP <- unlist(rp)
  }
  ##deal with empty RCLASS entries:
  else if(length(queried_data$RCLASS) == 0){
    queried_data$COMPOUNDS <- NA
    queried_data$RP <- NA
  }
  return(queried_data)
}
