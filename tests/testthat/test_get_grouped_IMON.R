library(PanViz)


testthat::test_that("tests that get_grouped_IMON produces an igraph object and throws correct errors on bad inputs", {
  df <- PanViz::GWAS_catalog_tsv_to_dataframe(file = system.file("extdata", "gwas-association-downloaded_2021-09-13-EFO_1000649.tsv", package="PanViz"))
  expect_equal(class(PanViz::get_grouped_IMON(dataframe = df, groupby = "studies", ego = 5, save_file = FALSE, colour_groups = FALSE)), "igraph")
  expect_equal(class(PanViz::get_grouped_IMON(dataframe = df, groupby = "studies", ego = 5, save_file = FALSE, colour_groups = TRUE)), "igraph")
  df <- data.frame(snps = c("rs3454594054594", "rs4545490954"), studies = c("study1", "study2"), traits = c("trait1", "trait1"))
  expect_error(object = PanViz::get_grouped_IMON(dataframe = df, groupby = "studies", ego = 5, save_file = FALSE), regexp =  "None of the supplied SNPs could be queried via dbSNP")
  df <- data.frame(snps = PanViz::er_snp_vector[1:2], studies = c("study1", "study2"), traits = c("trait1", "trait1"))
  expect_error(object = PanViz::get_grouped_IMON(dataframe = df, groupby = "studies", ego = 5, save_file = FALSE), regexp =  "None of the provided SNPs could be mapped to genes provided by the Kyoto Encyclopedia of Genes and Genomes")
  df <- PanViz::GWAS_catalog_tsv_to_dataframe(file = system.file("extdata", "gwas-association-downloaded_2021-09-13-EFO_1000649.tsv", package="PanViz"))
})
