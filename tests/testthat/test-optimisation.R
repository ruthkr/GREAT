# all_data_df <- system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "greatR") %>%
#   utils::read.csv() %>%
#   dplyr::filter(locus_name %in% c("BRAA05G005370.3C", "BRAA02G043220.3C"))
#
# # Test optimisation functions ----
#
# test_that("optimise_registration_params works", {
#   set.seed(123)
#   optimised_parameters <- optimise_registration_params(
#     input_df = all_data_df,
#     # genes = c("BRAA05G005370.3C", "BRAA02G043220.3C"),
#     initial_rescale = FALSE,
#     do_rescale = TRUE,
#     min_num_overlapping_points = 4,
#     accession_data_to_transform = "Col0",
#     accession_data_ref = "Ro18",
#     start_timepoint = "reference",
#     expression_value_threshold = 5,
#     is_data_normalised = FALSE,
#     num_iterations = 5,
#     boundary_coverage = 0.25
#   )
# })
#
# # registration_results_opt <- scale_and_register_data(
# #   input_df = data_test,
# #   stretches = 2.037,
# #   shifts = 0.62,
# #   min_num_overlapping_points = 4,
# #   initial_rescale = FALSE,
# #   do_rescale = TRUE,
# #   accession_data_to_transform = "Col0",
# #   accession_data_ref = "Ro18",
# #   start_timepoint = "reference",
# #   optimise_shift_extreme = FALSE
# # )
#
# # Preprocess data ----
#
# processed_data <- preprocess_data(
#   input_df = all_data_df,
#   initial_rescale = initial_rescale,
#   accession_data_to_transform = accession_data_to_transform,
#   accession_data_ref = accession_data_ref,
#   start_timepoint = start_timepoint,
#   expression_value_threshold = 5,
#   is_data_normalised = FALSE
# )
#
# all_data_df <- processed_data$all_data_df
# mean_df <- processed_data$mean_df
# mean_df_sc <- processed_data$mean_df_sc
# to_shift_df <- processed_data$to_shift_df
# time_to_add <- processed_data$time_to_add
#
# best_registration_list <- get_best_stretch_and_shift_after_optimisation(
#   to_shift_df,
#   all_data_df,
#   optimised_parameters = optimised_parameters,
#   do_rescale,
#   min_num_overlapping_points,
#   accession_data_to_transform,
#   accession_data_ref,
#   time_to_add
# )
#
# all_shifts <- best_registration_list$all_shifts
# best_shifts <- best_registration_list$best_shifts
# model_comparison_dt <- best_registration_list$model_comparison_dt
#
# shifted_mean_df <- apply_shift_to_registered_genes_only(
#   to_shift_df,
#   best_shifts,
#   model_comparison_dt,
#   accession_data_to_transform,
#   accession_data_ref,
#   time_to_add
# )
#
# imputed_mean_df <- impute_transformed_exp_values(
#   shifted_mean_df,
#   accession_data_to_transform,
#   accession_data_ref
# )
#
# imputed_mean_df %>%
#   plot_registration_results()
