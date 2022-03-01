#' Get BIC score from registering data - as part of optimisation
#'
#' @return
#' @export
get_BIC_from_registering_data <- function(input_df,
                                          stretches,
                                          shifts,
                                          min_num_overlapping_points,
                                          initial_rescale,
                                          do_rescale,
                                          accession_data_to_transform,
                                          accession_data_ref,
                                          start_timepoint = c("reference", "transform", "zero"),
                                          expression_value_threshold = 5,
                                          is_data_normalised = FALSE,
                                          optimise_shift_extreme = FALSE) {
  # Validate parameters
  start_timepoint <- match.arg(start_timepoint)

  # # Suppress "no visible binding for global variable" note
  # accession <- NULL
  # timepoint <- NULL
  # sc.expression_value <- NULL
  # expression_value <- NULL
  # locus_name <- NULL

  # Make sure the data are data.tables
  all_data_df <- data.table::as.data.table(input_df)

  mean_df <- get_mean_data(
    exp = all_data_df,
    expression_value_threshold = expression_value_threshold,
    accession_data_to_transform = accession_data_to_transform,
    is_data_normalised = is_data_normalised
  )

  # Parse start_timepoint
  if (start_timepoint == "reference") {
    time_to_add <- min(all_data_df[accession == accession_data_ref, timepoint])
  } else if (start_timepoint == "transform") {
    time_to_add <- min(all_data_df[accession == accession_data_to_transform, timepoint])
  } else {
    time_to_add <- 0
  }

  # Filter genes of original input data as in mean_df
  all_data_df <- all_data_df[all_data_df$locus_name %in% unique(mean_df$locus_name)]
  all_data_df <- subset(all_data_df, select = c("locus_name", "accession", "tissue", "timepoint", "expression_value", "group"))

  # TODO: validate colnames
  # Apply normalisation of expression for each gene across all timepoints
  mean_df_sc <- data.table::copy(mean_df)

  # Apply scaling
  mean_df_sc[, sc.expression_value := scale(expression_value, scale = TRUE, center = TRUE), by = .(locus_name, accession)]

  # Apply scaling before registration (if initial_rescale == TRUE), otherwise using original data
  if (initial_rescale == TRUE) {

    # apply rescale to mean_df prior to registration
    to_shift_df <- data.table::copy(mean_df_sc)
    to_shift_df$expression_value <- to_shift_df$sc.expression_value
    to_shift_df$sc.expression_value <- NULL

    # apply THE SAME rescale to all_data_df prior to registration
    all_data_df <- scale_all_rep_data(mean_df, all_data_df, "scale")
  } else {
    to_shift_df <- data.table::copy(mean_df)
  }

  L <- get_best_stretch_and_shift_simplified(
    to_shift_df,
    all_data_df,
    stretches,
    shifts,
    do_rescale,
    min_num_overlapping_points,
    accession_data_to_transform,
    accession_data_ref,
    time_to_add,
    optimise_shift_extreme
  )

  return(L$model_comparison_dt$registered.BIC - L$model_comparison_dt$separate.BIC.after)
}


#' Calculate best shifts and stretches for each gene, also calculate AIC/BIC under registration or non-registration (simplified)
#'
#' @noRd
get_best_stretch_and_shift_simplified <- function(to_shift_df,
                                                  all_data_df,
                                                  stretches,
                                                  shifts,
                                                  do_rescale,
                                                  min_num_overlapping_points,
                                                  accession_data_to_transform,
                                                  accession_data_ref,
                                                  time_to_add,
                                                  optimise_shift_extreme) {
  # # Suppress "no visible binding for global variable" note
  # is_best <- NULL
  # gene <- NULL
  # delta.BIC <- NULL
  #
  # Warning to make sure users have correct accession data
  if (!(accession_data_to_transform %in% all_data_df$accession & accession_data_ref %in% all_data_df$accession)) {
    stop("get_best_stretch_and_shift(): data accessions should have been converted to correct accession.")
  }


  # Calculate all the shift scores given this stretch. Score is mean(dist^2), over overlapping points if do_rescale=T, is rescaled by the mean FOR THE OVERLAPPING POINTS. (but not by the SD.)
  all_shifts <- calculate_all_best_shifts(
    mean_df = to_shift_df,
    stretch_factor = stretches,
    shifts,
    do_rescale,
    min_num_overlapping_points,
    accession_data_to_transform,
    accession_data_ref,
    optimise_shift_extreme
  )

  # Ensure no duplicated rows
  all_shifts <- unique(all_shifts)


  # Calculate the BIC & AIC for the best shifts found with this stretch.compared to treating the gene's expression separately in data to transform and reference data
  model_comparison_dt <- calculate_all_model_comparison_stats(
    all_data_df,
    all_shifts,
    accession_data_to_transform,
    accession_data_ref,
    time_to_add
  )

  return(
    list(model_comparison_dt = model_comparison_dt)
  )
}
