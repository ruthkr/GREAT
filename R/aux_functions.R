#' Get BIC score from registering data
#'
#' Simplified version of \code{\link{scale_and_register_data}} for \code{\link{optimise_registration_params}}.
#'
#' @noRd
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

  # Preprocess data
  processed_data <- preprocess_data(
    input_df = input_df,
    initial_rescale = initial_rescale,
    accession_data_to_transform = accession_data_to_transform,
    accession_data_ref = accession_data_ref,
    start_timepoint = start_timepoint,
    expression_value_threshold = expression_value_threshold,
    is_data_normalised = is_data_normalised
  )

  all_data_df <- processed_data$all_data_df
  to_shift_df <- processed_data$to_shift_df
  time_to_add <- processed_data$time_to_add

  # Calculate the best registration
  best_registration_list <- get_best_stretch_and_shift_simplified(
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

  registered_BIC <- best_registration_list$model_comparison_dt$registered.BIC
  separate_BIC <- best_registration_list$model_comparison_dt$separate.BIC
  BIC_diff <- registered_BIC - separate_BIC

  return(BIC_diff)
}

#' Calculate best shifts and stretches for each gene, also calculate BIC under registration or non-registration
#'
#' Simplified version of \code{\link{get_best_stretch_and_shift}} for \code{\link{optimise_registration_params}}.
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


  # Calculate the BIC for the best shifts found with this stretch.compared to treating the gene's expression separately in data to transform and reference data
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
