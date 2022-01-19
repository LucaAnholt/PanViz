#' GWAS Catalog tsv to cleaned dataframe
#'
#' @param file - The file name of a GWAS Catalog association .tsv file
#'
#' @return A processed dataframe containing only the columns including GWAS studies, traits/phenotypes and relevant SNPs in NCBI standard accession number naming convention
#' @export
#'
GWAS_catalog_tsv_to_dataframe <- function(file){
  if(missing(file)){ ##check if user has provided file, if not send error
    stop("Argument file missing. Please provide a directory path to a GWAS Catalog association .tsv file")
  }
  if(length(grep(pattern = ".tsv", x = file)) != 1){ ##check if user has provided correct file extension, if not send error
    stop("File provided must be .tsv file format")
  }
  file <- as.data.frame(data.table::fread(file)) #read tsv as dataframe
  if(is.null(file)){ ##check that file was properly read, if not send error
    stop(".tsv file is empty - please check the file provided")
  }
  results <- file %>%
    tidyr::separate_rows(data = ., SNPS, sep = ";\\s+") %>% #separate aggregated snps into separate rows
    tidyr::separate_rows(data = ., SNPS, sep = " x ") %>% #separate aggregated snps into separate rows
    ##only take snps in standard dbSNP accession number naming convention
    dplyr::filter(.data = ., stringr::str_detect(.data$SNPS, "rs")) %>%
    ##select SNPs, STUDY and DISEASE/TRAIT columns:
    dplyr::select(.data = ., SNPS, STUDY, `DISEASE/TRAIT`) %>%
    ##Rename these columns:
    dplyr::rename(.data = ., snps = SNPS, studies = STUDY, traits = `DISEASE/TRAIT`)
  return(results)
}
