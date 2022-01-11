id_table <- system.file("extdata/sample_data/id_table.rds", package = "greatR") %>% readRDS()
arabidopsis_expression <- system.file("extdata/sample_data/arabidopsis_expression.rds", package = "greatR") %>% readRDS()
brassica_rapa_expression <- system.file("extdata/sample_data/brassica_rapa_expression.rds", package = "greatR") %>% readRDS()

all_data_df <- system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "greatR") %>%
  utils::read.csv() %>%
  as.data.table()

mean_df <- system.file("extdata/brapa_arabidopsis_mean.csv", package = "greatR") %>%
  utils::read.csv() %>%
  as.data.table()

test_that("get_expression_of_interest works", {
  exp <- get_expression_of_interest(
    data_ref = brassica_rapa_expression,
    data_to_transform = arabidopsis_expression,
    id_table = id_table,
    lookup_col_ref_and_id_table = "CDS.model",
    lookup_col_ref_and_to_transform = "locus_name",
    colnames_wanted = NULL,
    tissue_wanted = "apex",
    gene_of_interest_acc = c("AT1G69120", "AT5G618"),
    sum_exp_data_ref = FALSE,
    accession_data_to_transform = "Col0"
  )

  expect_equal(class(exp)[[1]], "data.table")
})

test_that("get_mean_data works", {
  mean_data <- get_mean_data(
    all_data_df,
    expression_value_threshold = 5,
    accession_data_to_transform = "Col0"
  )

  expect_equal(
    round(mean_data$expression_value, 5),
    round(mean_df$expression_value, 5)
  )
})
