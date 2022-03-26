#' snp_gene_chr_match
#'
#' @param gene_loc - dataframe of genes and their chromosome numbers and start/stop positions
#' @param snp_loc - snp locations
#' @return - a recursive list of gene with their relative snps that have the same chromosome number
snp_gene_chr_match <- function(snp_loc, gene_loc){
  chr_match <- c(which(gene_loc$chr_n %in% snp_loc)) #find index of snps that have same chr number as input gene
  if(length(chr_match) != 0){ #checking gene format is correct
    return(chr_match)
  }
  else{
    return(NA)
  }
}
