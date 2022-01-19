library(PanViz)

testthat::test_that("tests that GWAS_catalog_tsv_to_dataframe produces a dataframe object and throws correct errors on bad inputs", {
  all(inherits(PanViz::GWAS_catalog_tsv_to_dataframe(file = paste0(find.package("PanViz", lib.loc=NULL, quiet = TRUE), "/data/gwas-association-downloaded_2021-09-13-EFO_1000649.tsv")), c("tbl_df", "data.frame"), which = TRUE) > 0)
  testthat::expect_error(PanViz::GWAS_catalog_tsv_to_dataframe(), regexp = "Argument file missing. Please provide a directory path to a GWAS Catalog association .tsv file")
  testthat::expect_error(PanViz::GWAS_catalog_tsv_to_dataframe(file = paste0(find.package("PanViz", lib.loc=NULL, quiet = TRUE), "/data/gwas-association-downloaded_2021-09-13-EFO_1000649.csv")), regexp = "File provided must be .tsv file format")
})
