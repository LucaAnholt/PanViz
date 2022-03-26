#' gene_cleanup
#'
#' @description Internal function for cleaning up queried KEGG gene data
#'
#' @param queried_data queried KEGG gene data (recursive list as returned from KEGGREST)
#'
#' @return cleaned up queried KEGG gene recursive lists, filtered by genes with adjacent enzymes
#'
gene_cleanup <- function(queried_data){
    ##deleting unnecessary data
    queried_data$DEFINITION <- NULL
    queried_data$ORGANISM <- NULL
    queried_data$PATHWAY <- NULL
    queried_data$BRITE <- NULL
    queried_data$MOTIF <- NULL
    queried_data$STRUCTURE <- NULL
    queried_data$AASEQ <- NULL
    queried_data$NTSEQ <- NULL
    queried_data$NETWORK <- NULL
    ##checking if gene has adjacent enzymes:
    if(all(stringr::str_detect(queried_data$ORTHOLOGY, "EC:")) & !is.null(queried_data$ORTHOLOGY)){
      enzyme <- stringr::str_extract_all(queried_data$ORTHOLOGY[[1]], "\\[EC:.*\\]")
      enzyme <- gsub("\\[", "", enzyme)
      enzyme <- gsub("\\]", "", enzyme)
      queried_data$enzyme <- enzyme
    }
    ##if no adjacent enzymes delete from recursive list
    else{
      queried_data <- NA
    }
    return(queried_data)
  }
