library(PanViz)

testthat::test_that("tests that get_IMON produces an igraph object and throws correct errors on bad inputs", {
  data("er_snp_vector")
  testthat::expect_equal(class(PanViz::get_IMON(snp_list = er_snp_vector, ego = 5, save_file = FALSE)), "igraph")
  testthat::expect_error(object = PanViz::get_IMON(snp_list = er_snp_vector[1:2], ego = 5, save_file = FALSE),regexp =  "None of the provided SNPs could be mapped to genes provided by the Kyoto Encyclopedia of Genes and Genomes")
  testthat::expect_error(object = PanViz::get_IMON(snp_list = c("rs230454923045"), ego = 5, save_file = FALSE), regexp =  "None of the supplied SNPs could be queried via dbSNP")
})

