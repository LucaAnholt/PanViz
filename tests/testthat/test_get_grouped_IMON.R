library(PanViz)


testthat::test_that("tests that get_grouped_IMON produces an igraph object and throws correct errors on bad inputs", {
  df <- PanViz::GWAS_data_reader(file = system.file("extdata", "gwas-association-downloaded_2021-09-13-EFO_1000649.tsv", package="PanViz"), snp_col = "SNPS", study_col = "STUDY", trait_col = "DISEASE/TRAIT")
  expect_equal(class(PanViz::get_grouped_IMON(dataframe = df, groupby = "studies", ego = 5, save_file = FALSE, colour_groups = FALSE)), "igraph")
  expect_equal(class(PanViz::get_grouped_IMON(dataframe = df, groupby = "studies", ego = 5, save_file = FALSE, colour_groups = TRUE)), "igraph")
  df <- data.frame(snps = c("rs3454594054594", "rs4545490954"), studies = c("study1", "study2"), traits = c("trait1", "trait1"))
  expect_error(object = PanViz::get_grouped_IMON(dataframe = df, groupby = "studies", ego = 5, save_file = FALSE), regexp =  "None of the supplied SNPs could be queried via dbSNP")
  data("er_snp_vector")
  df <- data.frame(snps = er_snp_vector[1:2], studies = c("study1", "study2"), traits = c("trait1", "trait1"))
  expect_error(object = PanViz::get_grouped_IMON(dataframe = df, groupby = "studies", ego = 5, save_file = FALSE), regexp =  "None of the provided SNPs could be mapped to genes provided by the Kyoto Encyclopedia of Genes and Genomes")
})
