#' Calculate extreme shifts after the stretch transformation
#'
#' `get_extreme_shifts_for_all` is used to calculate the minimum and maximum shifts can apply to data after the stretch transformation.
#'
#' @param mean_df Input dataframe containing data to transform and data fix.
#' @param stretch_factor Stretch transformation factor wanted.
#' @param min_num_overlapping_points Bound the extreme allowed shifts, such than at least this many timepoints are being compared for both accessions.
#' @param shift_extreme Approximation of maximum and minimum shifts allowed.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_fix Accession name of data fix.
#'
#' @return minimum and maximum values of shifts
#' @export
get_extreme_shifts_for_all <- function(mean_df,
                                       stretch_factor,
                                       min_num_overlapping_points,
                                       shift_extreme,
                                       accession_data_to_transform,
                                       accession_data_fix) {


  # This function is a wrapper for calc_extreme_shifts() to be able to move it out of the loop so don't calculate for every gene.
  # Cut dataframe to a single gene.
  curr_sym <- unique(mean_df$locus_name)[1]
  mean_df <- mean_df[mean_df$locus_name == curr_sym, ]

  # Transform timepoint to be time from first timepoint
  mean_df[, delta_time := timepoint - min(timepoint), by = .(accession)]

  # Apply stretch_factor to the data to transform, leave the data fix as it is
  mean_df$delta_time[mean_df$accession == accession_data_to_transform] <- mean_df$delta_time[mean_df$accession == accession_data_to_transform] * stretch_factor

  # Calculate max and min shifts
  max_min_shifts <- calc_extreme_shifts(
    mean_df,
    min_num_overlapping_points,
    shift_extreme,
    accession_data_to_transform,
    accession_data_fix
  )

  return(max_min_shifts)
}

#' Calculate extreme shifts
#'
#' `calc_extreme_shifts` is used to calculate the minimum and maximum shifts, whilst preserving the criteria that at least min_num_overlapping_points are being compared from both accessions.
#'
#' @param mean_df Input dataframe containing data to transform and data fix.
#' @param min_num_overlapping_points Bound the extreme allowed shifts, such than at least this many timepoints are being compared for both accessions.
#' @param shift_extreme Approximation of maximum and minimum shifts allowed.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_fix Accession name of data fix.
#'
#' @return minimum and maximum values of shifts
#'
#' @export
calc_extreme_shifts <- function(mean_df,
                                min_num_overlapping_points,
                                shift_extreme,
                                accession_data_to_transform = "Col0",
                                accession_data_fix = "Ro18") {

  # Make copy of the mean_df to make sure that the original dataframe will not edited
  original <- data.table::copy(mean_df)
  original$shifted_time <- original$delta_time


  # Negative extreme shift will be -1 * exactly the difference between 1 of the stretched data_to_transform (e.g. Col0) timepoints, and the smallest data_fix (e.g. R018) timepoint
  # Potitive extreme will be the difference between 1 of the ata_to_transform timepoints, and the maximum data_fix timepoint
  neg_extreme_candidate <- -1 * (original$delta_time[original$accession == accession_data_to_transform] - min(original$delta_time[original$accession == accession_data_fix]))
  pos_extreme_candidates <- max(original$delta_time[original$accession == accession_data_fix]) - original$delta_time[original$accession == accession_data_to_transform]

  # Among of these candidates, find the most extreme values which maintaining the required number of overlapping time-points to be considered.
  num_overlapping_points <- sapply(neg_extreme_candidate,
    FUN = calc_num_overlapping_points,
    data = original,
    accession_data_to_transform = accession_data_to_transform
  )
  if (all(num_overlapping_points < min_num_overlapping_points)) {
    stop(paste0(
      "calc_extreme_shifts():\nafter applying stretch factor:", stretch, " to ", transformed.timecourse, ", none of the considered shifts have ",
      "min_num_overlapping_points (", min_num_overlapping_points, ") overlapping timepoints with the other timecourse!\n",
      "maybe try a smaller stretch, and double check you're applying it to the correct timecourse."
    ))
  }
  neg_extreme <- min(neg_extreme_candidate[num_overlapping_points >= min_num_overlapping_points])

  num_overlapping_points <- sapply(pos_extreme_candidates,
                                   FUN = calc_num_overlapping_points,
                                   data = original,
                                   accession_data_to_transform = accession_data_to_transform)
  pos_extreme <- max(pos_extreme_candidates[num_overlapping_points >= min_num_overlapping_points])

  # Hard code maximum and minimum allowed shifts, as noticed spurious registrations when too extreme shifts allowed
  if (neg_extreme < (-1 * shift_extreme)) {
    neg_extreme <- -1 * shift_extreme
  }
  if (pos_extreme > 1 * shift_extreme) {
    pos_extreme <- shift_extreme
  }

  return(list(neg_extreme, pos_extreme))
}
