#' Full summary-level GWAS data for estrogen-receptor positive breast cancer
#' (meSH/HPO: D001943) downloaded from the GWAS Central on 8/2/2022
#'
#' A dataset containing a table of SNPs (summary-level GWAS data) associated
#' with estrogen-receptor positive breast cancer (meSH/HPO: D001943), and various other
#' study information, collated by the GWAS Central database.
#'
#' This data was sourced via the GWAS Central GWAS Mart API
#' (https://mart.gwascentral.org/), where the following filters and attributes
#' were applied:
#'
#' - Under "Filter", a filter for meSH/HPO Phenotype Identifier D001943 was applied,
#' corresponding to estrogen positive breast cancer.
#' - Under "Attributes", "Study Name" and "Annotation Name" was selected in STUDY
#' INFORMATION and HGNC Gene Symbol selected under ASSOCIATION RESULTS
#'
#' @docType data
#' @keywords external dataset
#' @references https://mart.gwascentral.org/
#' @usage PanViz::GWAS_data_reader(file = system.file("extdata", "GWASCentralMart_ERplusBC.tsv", package="PanViz"), snp_col = "Source Marker Accession", study_col = "Study Name", trait_col = "Annotation Name")
#' @format A tab separated file with 95 rows and 7 columns
"GWASCentralMart_ERplusBC.tsv"
