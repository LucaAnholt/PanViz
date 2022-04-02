#' get_compound_hashmap
#' @description This function creates a hashmap of associated KEGG IDs and compound/metabolite names
#' @return hashmap/recursive list of associated KEGG IDs and compound/metabolite names
get_compound_hashmap <- function(){
  compounds_Raw_IDs <- KEGGREST::keggList("compound")
  ##cleaning up raw data:
  Clean_Compounds <- attr(compounds_Raw_IDs, "names")
  ##only query compounds that don't exist in the database:
  Clean_Compounds <- gsub(pattern = "cpd:", replacement = "", x = Clean_Compounds)
  index <- which(Clean_Compounds %in% unique(unlist(unname(adjl_RP_C))))
  compound_names <- unname(compounds_Raw_IDs)
  compound_names <- lapply(X = compound_names, FUN = function(x){return(stringr::str_split(string = x, pattern = ";")[[1]][[1]])})
  names(compound_names) <- Clean_Compounds
  return(compound_names)
}
