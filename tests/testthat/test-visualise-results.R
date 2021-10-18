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
  accession_data_to_transform = "Col0",
  accession_data_ref = "Ro18",
  data_to_transform_time_added = 11,
  data_ref_time_added = 11
)


test_that("Function to plot of registration result returns gg/ggplot object", {
  class_plot_registered_gene_of_interest <- registration_results$imputed_mean_df %>%
    plot_registered_gene_of_interest() %>%
    class()

  expect_equal(class_plot_registered_gene_of_interest[1], "gg")
  expect_equal(class_plot_registered_gene_of_interest[2], "ggplot")
})

sample_distance_results <- calculate_between_sample_distance(
  registration_results$mean_df,
  registration_results$mean_df_sc,
  registration_results$imputed_mean_df,
  gene_col = "locus_name",
  compare_ref_vs_transform = TRUE,
  accession_data_ref = "Ro18"
)

test_that("Function make_heatmap to plot of distance heatmap of samples returns gg/ggplot object", {
  class_make_heatmap <- sample_distance_results$D.registered %>%
    make_heatmap(ylabel = "Distance") %>%
    class()

  expect_equal(class_make_heatmap[1], "gg")
  expect_equal(class_make_heatmap[2], "ggplot")
})
