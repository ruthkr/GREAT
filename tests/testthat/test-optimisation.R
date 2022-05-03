all_data_df <- system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "greatR") %>%
  utils::read.csv() %>%
  dplyr::filter(locus_name %in% c("BRAA05G005370.3C"))

# Pipeline parameters ----

stretches_bound <- c(2.0, 2.05)
shifts_bound <- c(0.60, 0.65)
initial_rescale <- FALSE
do_rescale <- TRUE
min_num_overlapping_points <- 4
accession_data_to_transform <- "Col0"
accession_data_ref <- "Ro18"
start_timepoint <- "reference"
expression_value_threshold <- 5
is_data_normalised <- FALSE

# Test optimisation functions ----

test_that("optimise_registration_params works", {
  set.seed(123)
  num_iterations <- 1

  optimised_parameters <- optimise_registration_params(
    input_df = all_data_df,
    stretches = stretches_bound,
    shifts = shifts_bound,
    initial_rescale = initial_rescale,
    do_rescale = do_rescale,
    min_num_overlapping_points = min_num_overlapping_points,
    accession_data_to_transform = accession_data_to_transform,
    accession_data_ref = accession_data_ref,
    start_timepoint = start_timepoint,
    expression_value_threshold = expression_value_threshold,
    is_data_normalised = is_data_normalised,
    num_iterations = num_iterations
  )

  processed_data <- preprocess_data(
    input_df = all_data_df,
    initial_rescale = initial_rescale,
    accession_data_to_transform = accession_data_to_transform,
    accession_data_ref = accession_data_ref,
    start_timepoint = start_timepoint,
    expression_value_threshold = expression_value_threshold,
    is_data_normalised = is_data_normalised
  )

  best_registration_list <- get_best_stretch_and_shift_after_optimisation(
    processed_data$to_shift_df,
    processed_data$all_data_df,
    optimised_parameters,
    do_rescale,
    min_num_overlapping_points,
    accession_data_to_transform,
    accession_data_ref,
    processed_data$time_to_add
  )

  # Expected output for optimise_registration_params()
  expected_optimum_params_df <- data.frame(
    stringsAsFactors = FALSE,
    gene = c("BRAA05G005370.3C"),
    stretch = c(2.045),
    shift = c(0.645),
    BIC_diff = c(-14.0280742834283),
    is_registered = c(TRUE)
  )
  expect_equal(names(optimised_parameters), c("optimum_params_df", "candidate_params_df"))
  expect_equal(optimised_parameters$optimum_params_df, expected_optimum_params_df)
  expect_true(optimised_parameters$optimum_params_df$is_registered)
  expect_equal(nrow(optimised_parameters$candidate_params_df), num_iterations)

  # Expected output for get_best_stretch_and_shift_after_optimisation()
  expect_equal(names(best_registration_list), c("all_shifts", "best_shifts", "model_comparison_dt"))
  expect_equal(ncol(best_registration_list$all_shifts), 8)
  expect_equal(ncol(best_registration_list$best_shifts), 9)
  expect_equal(ncol(best_registration_list$model_comparison_dt), 8)
})

test_that("get_boundary_box works", {
  boundary_box_manual <- get_boundary_box(
    input_df = all_data_df,
    stretches_bound = stretches_bound,
    shifts_bound = shifts_bound,
    accession_data_to_transform = accession_data_to_transform,
    accession_data_ref = accession_data_ref,
    min_num_overlapping_points = min_num_overlapping_points,
    expression_value_threshold = expression_value_threshold
  )

  boundary_box_auto <- get_boundary_box(
    input_df = all_data_df,
    stretches_bound = NA,
    shifts_bound = NA,
    accession_data_to_transform = accession_data_to_transform,
    accession_data_ref = accession_data_ref,
    min_num_overlapping_points = min_num_overlapping_points,
    expression_value_threshold = expression_value_threshold
  )

  # Expected output for get_boundary_box()
  expect_equal(round(unname(unlist(boundary_box_manual)), 2), c(2.02, 2, 2.05, 0.62, 0.6, 0.65))
  expect_equal(round(unname(unlist(boundary_box_auto)), 2), c(2.67, 1.3, 4.5, 0, -27, 19.8))
})
