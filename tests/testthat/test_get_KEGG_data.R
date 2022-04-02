library(PanViz)

testthat::test_that("tests that get_KEGG_data returns list object and throws correct errors on bad inputs", {
  testthat::expect_error(PanViz::get_Kegg_data(CPU = 5, sleep = 5, progress_bar = TRUE), regexp = "Can only use a maximum of 2 workers")
  ##get time sysdata.rda was last changed
  before_file_change <- file.info(system.file("R", "sysdata.rda", package="PanViz"))
  ##run function, that adds to sysdata.rda
  PanViz::get_Kegg_data(CPU = 1, sleep = 5, progress_bar = TRUE)
  ##see if file has successfully been altered
  after_file_change <- file.info(system.file("R", "sysdata.rda", package="PanViz"))
  boolean_result <- all(before_file_change == after_file_change)
  testthat::expect_false(object = boolean_result)
})
