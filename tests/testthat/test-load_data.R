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

test_that("Function to get all expression of interest dataframe/data table object", {
  id_table_path <- system.file("extdata/sample_data/id_table.RDS", package = "genereg")
  arabidopsis_expression_path <- system.file("extdata/sample_data/arabidopsis_expression.RDS", package = "genereg")
  brassica_rapa_expression_path <- system.file("extdata/sample_data/brassica_rapa_expression.RDS", package = "genereg")
  list_gene_of_interest <- c("AT1G69120", "AT5G61850", "AT2G45660")

  all_data <- genereg::get_expression_of_interest(brassica_rapa_expression_path,
                                                  arabidopsis_expression_path,
                                                  id_table_path,
                                                  tissue_wanted = 'apex',
                                                  curr_GoIs = list_gene_of_interest,
                                                  sum_brassicas = F) %>%
    class()

  expect_equal(all_data[1], "data.table")
  expect_equal(all_data[2], "data.frame")
})

test_that("Function to get mean df return a list of dataframe/data table object", {
  id_table_path <- system.file("extdata/sample_data/id_table.RDS", package = "genereg")
  arabidopsis_expression_path <- system.file("extdata/sample_data/arabidopsis_expression.RDS", package = "genereg")
  brassica_rapa_expression_path <- system.file("extdata/sample_data/brassica_rapa_expression.RDS", package = "genereg")
  list_gene_of_interest <- c("AT1G69120", "AT5G61850", "AT2G45660")

  mean_df <- genereg::load_mean_df(brassica_rapa_expression_path,
    arabidopsis_expression_path,
    id_table_path,
    tissue_wanted = "apex",
    curr_GoIs = list_gene_of_interest,
    sum_brassicas = F
  )

  expect_equal(mean_df %>% class, "list")
  expect_equal(mean_df[[1]] %>% class, "data.table")
  expect_equal(mean_df[[2]] %>% class, "data.frame")
})



