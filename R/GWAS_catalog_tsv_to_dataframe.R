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
  df <- tidyr::separate_rows(file, SNPS, sep = ";\\s+") #separate aggregated snps into separate rows
  df <- tidyr::separate_rows(df, SNPS, sep = " x ") #separate aggregated snps into separate rows
  snps <- df$SNPS[stringr::str_detect(df$SNPS, "rs")] #only take snps in standard dbSNP accession number naming convention
  studies <- df[stringr::str_detect(df$SNPS, "rs"),]$STUDY #get studies for these rows
  traits <- df[stringr::str_detect(df$SNPS, "rs"),]$`DISEASE/TRAIT` #get traits for these rows
  results <- data.frame(snps = snps, studies = studies, traits = traits) #create new dataframe
  results$snps <- as.character(results$snps) #convert data to strings isntead of factors
  results$studies <- as.character(results$studies)
  results$traits <- as.character(results$traits)
  return(results)
}
