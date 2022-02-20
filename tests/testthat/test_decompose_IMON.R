library(PanViz)

testthat::test_that("tests that decompose_IMON returns list object and throws correct errors on bad inputs", {
  data("er_snp_vector")
  G <- PanViz::get_IMON(snp_list = er_snp_vector)
  G_list <- PanViz::decompose_IMON(G = G)
  expect_equal(class(G_list), "list")
  expect_equal(length(G_list), 8)
  expect_equal(length(igraph::V(G_list[[1]])), 1780)
  expect_error(PanViz::decompose_IMON(), "Please provide an unnconnected graph (IMON)", fixed = TRUE)
  expect_error(PanViz::decompose_IMON(G_list[[1]]), "The IMON you provided is already fully connected! Cannot decompose!", fixed = TRUE)
})
