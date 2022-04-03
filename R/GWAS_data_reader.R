#' GWAS_data_reader
#'
#' @param file - Character (string) containing the directory path to a .tsv or
#' .csv file containing summary level GWAS data, typically this can be sourced
#' from major GWAS databases such as the GWAS Catalog or GWAS Central.
#'
#' @param snp_col - Character (string) reflecting the column name containing
#' the SNP (standard dbSNP accession number, e.g. rs992531) data. In data sourced
#' from the GWAS Catalog, this column will typically be named "SNPS" and in GWAS
#' Central this will typically be "Source Marker Accession".
#'
#' @param study_col - Character (string) reflecting the column name containing
#' the study names associated with each SNP. In data sourced from the GWAS Catalog,
#' this column will typically be named "STUDY" and in GWAS Central this will typically
#' be "Study Name".
#'
#' @param trait_col - Character (string) reflecting the column name containing
#' the trait/phenotype names associated with each SNP. In data sourced from the
#' GWAS Catalog, this column will typically be named "DISEASE/TRAIT" and in GWAS
#' Central this will typically be "Annotation Name".
#'
#' @return A processed dataframe containing only the columns including GWAS studies,
#' traits/phenotypes and relevant SNPs in NCBI standard accession number naming
#' convention
#'
#' @examples
#' ##getting directory path to GWAS Catalog association .tsv file:
#' path = system.file("extdata",
#'   "gwas-association-downloaded_2021-09-13-EFO_1000649.tsv",
#'   package="PanViz")
#' ##opening/cleaning data:
#' df <- PanViz::GWAS_data_reader(file = path,
#'   snp_col = "SNPS",
#'   study_col = "STUDY",
#'   trait_col = "DISEASE/TRAIT")
#' ##getting directory path to GWAS Central association .tsv file:
#' path = system.file("extdata", "GWASCentralMart_ERplusBC.tsv",
#'   package="PanViz")
#' ##opening/cleaning data:
#' df <- PanViz::GWAS_data_reader(file = path,
#'   snp_col = "Source Marker Accession",
#'   study_col = "Study Name",
#'   trait_col = "Annotation Name")
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
#'
GWAS_data_reader <- function(file, snp_col, study_col, trait_col){
  ##check if user has provided file, if not throw error
  if(missing(file)){
    stop("Argument file missing. Please provide a directory path to a .tsv file containg summary-level GWAS data")
  }
  ##check if user has provided correct file extension, if not throw error
  file_extensions <- c(".tsv", ".csv")
  if(!any(stringr::str_detect(string = file, pattern = file_extensions)) == TRUE){
    stop("File provided must be .tsv or .csv file format")
  }
  ##read tsv as dataframe
  file <- tibble::tibble(data.table::fread(file))
  ##check that file was properly read, if not throw error
  if(is.null(file)){
    stop("The file provided could not be read - please check the file provided")
  }
  ##check that provided column names are provided as characters, if provided as
  ##indexes, then convert to characters:
  if(is.numeric(snp_col)){
    snp_col <- colnames(file)[[snp_col]]
  }
  if(is.numeric(study_col)){
    study_col <- colnames(file)[[study_col]]
  }
  if(is.numeric(trait_col)){
    trait_col <- colnames(file)[[trait_col]]
  }
  ##separating SNP rows, if necessary:
  file <- tidyr::separate_rows(file, all_of(snp_col), convert = TRUE)
  ##check that standard dbSNP accession number naming convention exists within the snp column:
  ##first check if the standard acession naming exists at all:
  if(!any(stringr::str_detect(string = as.character(file[[snp_col]]), pattern = "rs")) == TRUE){
    stop("The file and SNP column identifier provided does not return the expected standard dbSNP accession number e.g. 'rs73370840'")
  }
  ##next check if it does exist, if there also exists unexpected non-standard naming too (if so filter only standard naming)
  if(!any(stringr::str_detect(string = as.character(file[[snp_col]]), pattern = "rs")) == FALSE & all(stringr::str_detect(string = as.character(file[[snp_col]]), pattern = "rs")) == FALSE){
    file <- dplyr::filter(file, stringr::str_detect(string = as.character(file[[snp_col]]), pattern = "rs"))
  }
  ##select and rename required columns:
  file <- file %>%
    dplyr::select(all_of(snp_col), all_of(study_col), all_of(trait_col)) %>%
    dplyr::rename(snps = snp_col, studies = study_col, traits = trait_col)
  return(data.frame(file))
}
