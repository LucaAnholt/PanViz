#' NCBI Gene Location
#'
#' @param Gene_Enterez_IDs NCBI (enterez) Gene Ids to query
#'
#' @return dataframe containing gene IDs, chromosome number and chromosome start/stop sites
#'
NCBI_Gene_Locations <- function(Gene_Enterez_IDs){
  split_data <- split(Gene_Enterez_IDs, ceiling(seq_along(Gene_Enterez_IDs)/100))
  #Getting gene locations:
  cat("Querying gene locations from NCBI\n")
  pb <- utils::txtProgressBar(max = length(split_data), style = 3)
  raw_data_list <- list()
  genomic_data <- list()
  for(i in 1:length(split_data)){
    raw_data_list[[i]] <- rentrez::entrez_summary(db = "gene", id = split_data[[i]])
    pctg <- paste(round(i/length(split_data) *100, 0), "% completed")
    utils::setTxtProgressBar(pb, i, label = pctg)
  }
  close(pb)
  ##removing recursion in list:
  raw_data <- unlist(raw_data_list, recursive = FALSE)
  ##removing any genes that are not listed in the NCBI database
  raw_data <- lapply(raw_data, NCBI_clean)
  raw_data <- raw_data[!vapply(X = raw_data, FUN = function(x) all(is.na(x)), FUN.VALUE = as.logical(1))]
  ##accessing the genomic location data:
  raw_data <- lapply(raw_data, NCBI_clean_2)
  ##organising genomic information as dataframe:
  gene_loc <- as.data.frame(do.call(rbind, raw_data))
  ##add in gene IDs:
  gene_loc$gene_id <- rownames(gene_loc)
  ##rename columns:
  colnames(gene_loc) <- c("chr_n", "chr_start", "chr_stop", "gene_id")
  ##remove any rows with NA:
  gene_loc <- stats::na.omit(gene_loc)
  return(gene_loc)
}


