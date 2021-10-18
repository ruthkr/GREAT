all_data_df <- system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "GREAT") %>%
  utils::read.csv()

mean_df <- system.file("extdata/brapa_arabidopsis_mean.csv", package = "GREAT") %>%
  utils::read.csv()


registration_results <- scale_and_register_data(
  mean_df = mean_df,
  all_data_df = all_data_df,
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



test_that("calculate_between_sample_distance works", {

  sample_distance_results <- calculate_between_sample_distance(registration_results$mean_df,
                                                                registration_results$mean_df_sc,
                                                                registration_results$imputed_mean_df)

  expect_equal(class(sample_distance_results), "list")
  expect_equal(length(names(sample_distance_results)), 6)
  expect_equal(class(sample_distance_results$D.mean)[[1]], "data.table")
  expect_equal(class(sample_distance_results$D.scaled)[[1]], "data.table")
  expect_equal(class(sample_distance_results$D.registered)[[1]], "data.table")
  expect_equal(class(sample_distance_results$D.scaled.onlyNR)[[1]], "data.table")
  expect_equal(class(sample_distance_results$D.scaled.onlyR)[[1]], "data.table")
  expect_equal(class(sample_distance_results$D.registered.onlyR)[[1]], "data.table")

})
