library(PanViz)

testthat::test_that("tests that decompose_IMON returns list object and throws correct errors on bad inputs", {
  path <- system.file("testdata", "testIMON.Rds",package="PanViz")
  if(path == "" | is.null(path)){
    path <- "inst/testdata/testIMON.Rds"
  }
  G <- readRDS(path)
  G_list <- PanViz::decompose_IMON(G = G)
  expect_equal(class(G_list), "list")
  expect_equal(length(G_list), 14)
  expect_equal(length(igraph::V(G_list[[1]])), 2775)
  expect_error(PanViz::decompose_IMON(), "Please provide an unnconnected graph (IMON)", fixed = TRUE)
  expect_error(PanViz::decompose_IMON(G_list[[1]]), "The IMON you provided is already fully connected! Cannot decompose!", fixed = TRUE)
})
