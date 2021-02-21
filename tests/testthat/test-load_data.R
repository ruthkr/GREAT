test_that("Function to combine all data return dataframe/data table object", {
  id_table_path <- system.file("extdata/sample_data/id_table.RDS", package = "genereg")
  arabidopsis_expression_path <- system.file("extdata/sample_data/arabidopsis_expression.RDS", package = "genereg")
  brassica_rapa_expression_path <- system.file("extdata/sample_data/brassica_rapa_expression.RDS", package = "genereg")

  all_data <- genereg::get_all_data(brassica_rapa_expression_path,
                                           arabidopsis_expression_path,
                                           id_table_path) %>%
    class()

  expect_equal(all_data[1], "data.table")
  expect_equal(all_data[2], "data.frame")
})
