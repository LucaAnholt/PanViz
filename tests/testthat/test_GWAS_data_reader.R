library(PanViz)

testthat::test_that("tests that GWAS_data_reader produces a dataframe object and throws correct errors on bad inputs", {
  test1 <- system.file("extdata", "gwas-association-downloaded_2021-09-13-EFO_1000649.tsv", package="PanViz")
  test2 <- system.file("extdata", "GWASCentralMart_ERplusBC.tsv", package="PanViz")
  expect_equal(class(PanViz::GWAS_data_reader(file = test1, snp_col = "SNPS", study_col = "STUDY", trait_col = "DISEASE/TRAIT")), "data.frame")
  expect_equal(dim(PanViz::GWAS_data_reader(file = test1, snp_col = "SNPS", study_col = "STUDY", trait_col = "DISEASE/TRAIT")), c(110, 3))
  expect_equal(colnames(PanViz::GWAS_data_reader(file = test1, snp_col = "SNPS", study_col = "STUDY", trait_col = "DISEASE/TRAIT")), c("snps", "studies", "traits"))
  expect_equal(class(PanViz::GWAS_data_reader(file = test2, snp_col = "Source Marker Accession", study_col = "Study Name", trait_col = "Annotation Name")), "data.frame")
  expect_equal(dim(PanViz::GWAS_data_reader(file = test2, snp_col = "Source Marker Accession", study_col = "Study Name", trait_col = "Annotation Name")), c(95, 3))
  expect_equal(colnames(PanViz::GWAS_data_reader(file = test2, snp_col = "Source Marker Accession", study_col = "Study Name", trait_col = "Annotation Name")), c("snps", "studies", "traits"))
  expect_error(PanViz::GWAS_data_reader(file = "test.file", snp_col = "test", study_col = "test", trait_col = "test"), "File provided must be .tsv or .csv file format", fixed = T)
})

