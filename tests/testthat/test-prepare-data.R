id_table_path <- system.file("extdata/sample_data/id_table.RDS", package = "GREAT")
arabidopsis_expression_path <- system.file("extdata/sample_data/arabidopsis_expression.RDS", package = "GREAT")
brassica_rapa_expression_path <- system.file("extdata/sample_data/brassica_rapa_expression.RDS", package = "GREAT")

test_that("get_mean_and_all_exp_data works", {
  test_load_all_and_mean_df <- get_mean_and_all_exp_data(
    filepath_data_ref = brassica_rapa_expression_path,
    filepath_data_to_transform = arabidopsis_expression_path,
    filepath_id_table = id_table_path,
    fix_id_table_shared_colname = "CDS.model",
    fix_and_to_transform_data_shared_colname = "locus_name",
    colnames_id_table = c("CDS.model", "symbol", "locus_name"),
    colnames_wanted = NULL,
    tissue_wanted = "apex",
    curr_GoIs = c("AT1G69120", "AT5G618"),
    sum_exp_data_ref = FALSE,
    accession_data_to_transform = "Col0",
    ids_data_ref_colnames = c("CDS.model", "locus_name"),
    max_expression_value_wanted = 5,
    exp_threshold = 0.5
  )

  expect_equal(class(test_load_all_and_mean_df), "list")
  # expect_equal(length(names(test_load_all_and_mean_df)), 2)
  expect_equal(class(test_load_all_and_mean_df[[1]])[[1]], "data.table")
  expect_equal(class(test_load_all_and_mean_df[[2]])[[1]], "data.table")
})
