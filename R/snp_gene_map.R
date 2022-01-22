#' Fast vectorised SNP to gene chromosome number and genomic location mapping
#'
#' @param gene_loc dataframe containing KEGG genes and relevant chromosome number and positions
#' @param snp_loc dataframe containing queried SNPs and relevant chromosome number and positions
#'
#' @return an adjacency list of SNPs with their relevant mapped genes to their genomic location
snp_gene_map <- function(gene_loc, snp_loc){
  ##finding genes and snps that match in chromosome number
  chr_matches <- lapply(snp_loc$chr_n, snp_gene_chr_match, gene_loc)
  ##finding matches in location for all chromosome matches (within a range of 1 million nucleotides):
  snp_gene_adj <- function(chr_matches, snp_loc, gene_loc, index){
    boolean_mask <- data.table::between(snp_loc[index, ]$chr_pos,
                            unname(unlist(gene_loc[chr_matches[[index]],]$chr_start)) - 500000,
                            unname(unlist(gene_loc[chr_matches[[index]],]$chr_stop)) + 500000)
    if(!any(boolean_mask)){
      return(NA)
    }
    else{
      mapped_genes <- gene_loc[chr_matches[[index]][boolean_mask], ]$gene_id
      return(mapped_genes)
    }
  }
  for(i in seq_along(chr_matches)){
    chr_matches[[i]] <- snp_gene_adj(chr_matches = chr_matches, snp_loc = snp_loc, gene_loc = gene_loc, index = i)
  }
  ##naming adjacency list:
  names(chr_matches) <- paste0("rs", as.character(snp_loc$snp_id))
  ##removing NA element values
  chr_matches <- chr_matches[!is.na(chr_matches)]
  return(chr_matches)
}
