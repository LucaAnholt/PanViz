library(PanViz)

testthat::test_that("tests that get_KEGG_data returns list object and throws correct errors on bad inputs", {
  testthat::expect_error(PanViz::get_Kegg_data(CPU = 5, sleep = 5), regexp = "Can only use a maximum of 2 workers")
  ##get time sysdata.rda was last changed
  before_file_change <- file.info("R/sysdata.rda")
  ##run function, that adds to sysdata.rda
  PanViz::get_Kegg_data(CPU = 2, sleep = 5)
  ##see if file has successfully been altered
  after_file_change <- file.info("R/sysdata.rda")
  testthat::expect_equal(object = all(before_file_change == after_file_change), expected = FALSE)
})
