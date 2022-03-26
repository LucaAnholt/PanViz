#' Full summary-level GWAS data for estrogen-receptor positive breast cancer
#' (EFO_1000649) downloaded from the GWAS Catalog on 13/9/2021
#'
#' A dataset containing a table of SNPs (summary-level GWAS data) associated
#' with estrogen-receptor positive breast cancer (EFO_1000649), and various other
#' study information, collated by the GWAS Catalog. This was downloaded via the
#' "Download Catalog data" link at https://www.ebi.ac.uk/gwas/efotraits/EFO_1000649.
#'
#' @docType data
#' @keywords external dataset
#' @references https://www.ebi.ac.uk/gwas/efotraits/EFO_1000649
#' @usage PanViz::GWAS_data_reader(file = system.file("extdata", "gwas-association-downloaded_2021-09-13-EFO_1000649.tsv", package="PanViz"), snp_col = "SNPS", study_col = "STUDY", trait_col = "DISEASE/TRAIT")
#' @format A tab separated file with 112 rows and 38 columns
"gwas-association-downloaded_2021-09-13-EFO_1000649.tsv"
