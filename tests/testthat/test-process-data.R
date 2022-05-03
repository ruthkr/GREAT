all_data_df <- system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "greatR") %>%
  utils::read.csv()

# Test preliminary analysis functions ----

test_that("scale_and_register_data works", {
  est_stretch <- get_approximate_stretch(
    all_data_df,
    accession_data_to_transform = "Col0",
    accession_data_ref = "Ro18"
  )

  expect_equal(ceiling(est_stretch), 3)
  expect_equal(floor(est_stretch), 2)
})


# Test scale_and_register_data() ----

test_that("scale_and_register_data works", {
  # Check scale_and_register_data()
  registration_results <- scale_and_register_data(
    input_df = all_data_df,
    stretches = c(2, 1.5, 1),
    shifts = seq(-4, 4, length.out = 27),
    min_num_overlapping_points = 4,
    maintain_min_num_overlapping_points = TRUE,
    initial_rescale = FALSE,
    do_rescale = TRUE,
    accession_data_to_transform = "Col0",
    accession_data_ref = "Ro18",
    start_timepoint = "reference"
  )

  expect_equal(class(registration_results), "list")
  expect_equal(class(registration_results$mean_df)[[1]], "data.table")
  expect_equal(class(registration_results$mean_df_sc)[[1]], "data.table")
  expect_equal(class(registration_results$imputed_mean_df)[[1]], "data.table")
  expect_equal(class(registration_results$all_shifts_df)[[1]], "data.table")
  expect_equal(class(registration_results$model_comparison_df)[[1]], "data.table")

  expected_model_comparison_df <- data.table::data.table(
    gene = c("BRAA02G043220.3C", "BRAA03G023790.3C", "BRAA03G051930.3C", "BRAA05G005370.3C", "BRAA07G034100.3C", "BRAA09G045310.3C", "BRAA02G018970.3C", "BRAA04G005470.3C", "BRAA06G025360.3C", "BRAA07G030470.3C"),
    separate.BIC = c(127.91726, 140.71107, 175.19107, 136.12161, 70.77961, 82.38534, 75.75246, 75.92666, 129.71473, 76.36370),
    stretch = c(2, 2, 2, 2, 2, 2, 1.5, 1.5, 1.5, 1.5),
    shift = c(-3.69230, 0.61538, 3.07692, 0.61538, -4, 0, -0.30769, 4, 1.84615, -0.30769)
  )

  expect_equal(
    registration_results$model_comparison_df$gene,
    expected_model_comparison_df$gene
  )
  expect_equal(
    round(registration_results$model_comparison_df$separate.BIC, 2),
    round(expected_model_comparison_df$separate.BIC, 2)
  )
  expect_equal(
    round(registration_results$model_comparison_df$stretch, 2),
    round(expected_model_comparison_df$stretch, 2)
  )
  expect_equal(
    round(registration_results$model_comparison_df$shift, 2),
    round(expected_model_comparison_df$shift, 2)
  )

  # Check calculate_between_sample_distance()
  sample_distance_results <- calculate_between_sample_distance(
    registration_results,
    gene_col = "locus_name",
    compare_ref_vs_transform = TRUE,
    accession_data_ref = "Ro18"
  )

  expect_equal(class(sample_distance_results), "list")
  expect_equal(length(names(sample_distance_results)), 6)
  expect_equal(class(sample_distance_results$distance_mean_df)[[1]], "data.table")
  expect_equal(class(sample_distance_results$distance_scaled_mean_df)[[1]], "data.table")
  expect_equal(class(sample_distance_results$distance_registered_df)[[1]], "data.table")
  expect_equal(class(sample_distance_results$distance_scaled_mean_df_only_nonreg)[[1]], "data.table")
  expect_equal(class(sample_distance_results$distance_scaled_mean_df_only_reg)[[1]], "data.table")
  expect_equal(class(sample_distance_results$distance_registered_df_only_reg)[[1]], "data.table")

  # Check plots
  gg_registered <- registration_results$imputed_mean_df %>%
    plot_registration_results()

  expect_equal(class(gg_registered)[[1]], "gg")
  expect_equal(class(gg_registered)[[2]], "ggplot")

  gg_distance <- sample_distance_results$distance_registered_df %>%
    plot_heatmap()

  expect_equal(class(gg_distance)[[1]], "gg")
  expect_equal(class(gg_distance)[[2]], "ggplot")

  # Check summary
  summary <- summary_model_comparison(registration_results$model_comparison_df)

  expect_equal(length(summary), 3)
  expect_equal(class(summary$df_summary), "data.frame")
  expect_equal(class(summary$registered_genes), "character")
  expect_equal(class(summary$non_registered_genes), "character")
})


# Pipeline parameters ----

stretches <- c(2, 1.5, 1)
shifts <- seq(-4, 4, length.out = 27)
min_num_overlapping_points <- 4
initial_rescale <- FALSE
do_rescale <- TRUE
accession_data_to_transform <- "Col0"
accession_data_ref <- "Ro18"
time_to_add <- 11

# Test each function in function: scale_and_register_data() ----

test_that("all functions called in scale_and_register_data work", {
  # Make sure the data are data.tables
  # mean_df <- data.table::as.data.table(mean_df)
  all_data_df <- data.table::as.data.table(all_data_df)

  mean_df <- get_mean_data(
    exp = all_data_df,
    expression_value_threshold = 5,
    accession_data_to_transform = accession_data_to_transform
  )

  # Apply normalisation of expression for each gene across all timepoints
  mean_df_sc <- data.table::copy(mean_df)

  # Apply scaling
  mean_df_sc[, sc.expression_value := scale(expression_value, scale = TRUE, center = TRUE), by = .(locus_name, accession)]

  # Apply scaling before registration (if initial_rescale == TRUE), otherwise using original data
  # if (initial_rescale == TRUE)
  # Apply rescale to mean_df prior to registration
  to_shift_df <- data.table::copy(mean_df_sc)
  to_shift_df$expression_value <- to_shift_df$sc.expression_value
  to_shift_df$sc.expression_value <- NULL

  # apply THE SAME rescale to all_data_df prior to registration
  all_data_df_scaled <- scale_all_rep_data(mean_df, all_data_df, "scale")

  # Expected output for scale_all_rep_data()
  expect_equal(class(all_data_df_scaled)[[1]], "data.table")
  expect_equal(nrow(all_data_df_scaled), nrow(all_data_df))

  # else (initial_rescale == FALSE)
  to_shift_df <- data.table::copy(mean_df)

  # Calculate the best registration. Returns all tried registrations, best stretch and shift combo,
  # and BIC stats for comparison of best registration model to separate models for expression of
  # each gene in Ro18 and Col0
  best_registration_list <- get_best_stretch_and_shift(
    to_shift_df,
    all_data_df,
    stretches,
    shifts,
    do_rescale,
    min_num_overlapping_points,
    maintain_min_num_overlapping_points = TRUE,
    accession_data_to_transform,
    accession_data_ref,
    time_to_add
  )

  all_shifts <- best_registration_list$all_shifts
  best_shifts <- best_registration_list$best_shifts
  model_comparison_dt <- best_registration_list$model_comparison_dt

  # Expected output for get_best_stretch_and_shift()
  expect_equal(class(best_registration_list), "list")
  expect_equal(names(best_registration_list), c("all_shifts", "best_shifts", "model_comparison_dt"))
  expect_equal(class(all_shifts)[[1]], "data.table")
  expect_equal(class(best_shifts)[[1]], "data.table")
  expect_equal(class(model_comparison_dt)[[1]], "data.table")

  # Report model comparison results
  # Get the best-shifted and stretched mean gene expression, only to genes which registration is better than
  # separate models by BIC. Don't stretch out, or shift genes for which separate is better.
  shifted_mean_df <- apply_shift_to_registered_genes_only(
    to_shift_df,
    best_shifts,
    model_comparison_dt,
    accession_data_to_transform,
    accession_data_ref,
    time_to_add
  )

  # Expected output for apply_shift_to_registered_genes_only()
  expect_equal(class(shifted_mean_df)[[1]], "data.table")
  expect_equal(nrow(shifted_mean_df), nrow(mean_df))
  expect_equal(class(shifted_mean_df$is_registered), "logical")

  # Impute transformed values at times == to the observed reference data points for each shifted transformed gene so can compare using heat maps.
  # transformed curves are the ones that been shifted around. Linear impute values for these
  # curves so that reference data samples can be compared to an transformed data point.
  imputed_mean_df <- impute_transformed_exp_values(
    shifted_mean_df,
    accession_data_to_transform,
    accession_data_ref
  )

  # Expected output for apply_shift_to_registered_genes_only()
  expect_equal(class(imputed_mean_df)[[1]], "data.table")
  expect_equal(class(imputed_mean_df$is_registered), "logical")
  # expect_gte(nrow(imputed_mean_df), nrow(mean_df))
  # expect_lte(nrow(imputed_mean_df), nrow(all_data_df))
  expect_equal(colnames(imputed_mean_df), colnames(shifted_mean_df))

  out <- list(
    "mean_df" = mean_df,
    "mean_df_sc" = mean_df_sc,
    "imputed_mean_df" = imputed_mean_df,
    "all_shifts" = all_shifts,
    "model_comparison_dt" = model_comparison_dt
  )
})
