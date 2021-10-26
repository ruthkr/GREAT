id_table <- system.file("extdata/sample_data/id_table.RDS", package = "GREAT") %>% readRDS()
arabidopsis_expression <- system.file("extdata/sample_data/arabidopsis_expression.RDS", package = "GREAT") %>% readRDS()
brassica_rapa_expression <- system.file("extdata/sample_data/brassica_rapa_expression.RDS", package = "GREAT") %>% readRDS()

all_data_df <- system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "GREAT") %>%
  utils::read.csv() %>%
  as.data.table()
mean_df <- system.file("extdata/brapa_arabidopsis_mean.csv", package = "GREAT") %>%
  utils::read.csv() %>%
  as.data.table()


test_that("get_expression_of_interest works", {
exp <- get_expression_of_interest(
      data_ref = brassica_rapa_expression,
      data_to_transform = arabidopsis_expression,
      id_table = id_table,
      fix_id_table_shared_colname = "CDS.model",
      fix_and_to_transform_data_shared_colname = "locus_name",
      colnames_id_table = c("CDS.model", "symbol", "locus_name"),
      colnames_wanted = NULL,
      tissue_wanted = "apex",
      curr_GoIs = c("AT1G69120", "AT5G618"),
      sum_exp_data_ref = FALSE,
      accession_data_to_transform = "Col0",
      ids_data_ref_colnames = c("CDS.model", "locus_name")
      )

  expect_equal(class(exp)[[1]], "data.table")
})

test_that("get_mean_data works", {
  mean_data <- get_mean_data(all_data_df,
              max_expression_value_wanted = 5,
              accession_data_to_transform = "Col0")

  expect_equal(
    round(mean_data$expression_value, 5),
    round(mean_df$expression_value, 5)
  )
})
