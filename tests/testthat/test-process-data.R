data_all <- system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "GREAT") %>%
  utils::read.csv()

data_mean <- system.file("extdata/brapa_arabidopsis_mean.csv", package = "GREAT") %>%
  utils::read.csv()

test_that("scale_and_register_data works", {
  registration_results <- scale_and_register_data(
    mean_df = data_mean,
    all_data_df = data_all,
    stretches = c(2, 1.5, 1),
    shift_extreme = 4,
    num_shifts = 27,
    min_num_overlapping_points = 4,
    initial_rescale = FALSE,
    do_rescale = TRUE,
    testing = FALSE,
    accession_data_to_transform = "Col0",
    accession_data_ref = "Ro18",
    data_to_transform_time_added = 11,
    data_ref_time_added = 11
  )

  expect_equal(class(registration_results), "list")
  expect_equal(class(registration_results$mean_df)[[1]], "data.table")
  expect_equal(class(registration_results$mean_df_sc)[[1]], "data.table")
  expect_equal(class(registration_results$imputed_mean_df)[[1]], "data.table")
  expect_equal(class(registration_results$all_shifts)[[1]], "data.table")
  expect_equal(class(registration_results$model_comparison_dt)[[1]], "data.table")
})
