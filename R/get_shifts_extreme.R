# mean_df <- mean_df
# stretch_factor
# min_num_overlapping_points
# shift_extreme
#' @export
get_extreme_shifts_for_all <- function(mean_df, stretch_factor, min_num_overlapping_points, shift_extreme) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # wrapper for calc_extreme_shifts to be able to move it out of the loop so don't calculate for every gene.

  #min_num_overlapping_points <- 5 # bound the extreme allowed shifts, such than at least this many timepoints are being compared for both accessions.
  # cut data.table to a single gene
  curr_sym <- unique(mean_df$locus_name)[1]
  mean_df <- mean_df[mean_df$locus_name==curr_sym, ]

  # transform timepoint to be time from first timepoint
  mean_df[, delta_time:=timepoint - min(timepoint), by=.(accession)]
  # apply stretch_factor to the arabidopsis, leave the rapa as is
  mean_df$delta_time[mean_df$accession=='Col0'] <- mean_df$delta_time[mean_df$accession=='Col0']*stretch_factor

  # calculate min shift and max time shift, which still allows overlap of at least 5 times to be compared from whichever accession will be considering fewer timepoints from.
  # Shift is applied to the arabidopsis - so the 5th largest arabidopsis time is the biggest -ve shift can be applied
  # and the biggest shift which can be applied is to make the 5th smallest arabidopsis time == largest brassica time
  #data.table::setorder(mean_df, delta_time)
  #min_shift <- min(mean_df$delta_time[mean_df$accession=='Ro18']) - mean_df$delta_time[mean_df$accession=='Col0'][length(mean_df$delta_time[mean_df$accession=='Col0'])-4]
  #max_shift <- max(mean_df$delta_time[mean_df$accession=='Ro18']) - mean_df$delta_time[mean_df$accession=='Col0'][5]

  M <- calc_extreme_shifts(mean_df, min_num_overlapping_points, shift_extreme)
  return(M)
}

#' Calculate extreme shifts
#'
#' `calc_extreme_shifts` is used to calculate the minimum and maximum shifts can apply to data after the stretch transformation, whilst preserving the criteria that at least min_num_overlapping_points are being compared from both accessions.
#'
#' @param mean_df Input dataframe containing data to align and data target.
#' @param min_num_overlapping_points Minimum number of overlapping points allowed.
#' @param shift_extreme Approximation of maximum and minimum shifts allowed.
#' @param accession_data_to_align Accession name of data which will be aligned.
#' @param accession_data_target Accession name of data target.
#'
#' @return minimum and maximum values of shifts
#'
#' @export
calc_extreme_shifts <- function(mean_df,
                                min_num_overlapping_points,
                                shift_extreme,
                                accession_data_to_align = "Col0",
                                accession_data_target = "Ro18") {

  # Make copy of the mean_df to make sure that the original dataframe will not edited
  original <- data.table::copy(mean_df)
  original$shifted_time <- original$delta_time


  # Negative extreme shift will be -1 * exactly the difference between 1 of the stretched data_to_align (e.g. Col0) timepoints, and the smallest data_target (e.g. R018) timepoint
  # Potitive extreme will be the difference between 1 of the ata_to_align timepoints, and the maximum data_target timepoint
  neg_extreme_candidate <- -1 * (original$delta_time[original$accession == accession_data_to_align] - min(original$delta_time[original$accession == accession_data_target]))
  pos_extreme_candidates <- max(original$delta_time[original$accession == accession_data_target]) - original$delta_time[original$accession == accession_data_to_align]

  # Among of these candidates, find the most extreme values which maintaining the required number of overlapping time-points to be considered.
  num_overlapping_points <- sapply(neg_extreme_candidate,
    FUN = calc_num_overlapping_points,
    original = original
  )
  if (all(num_overlapping_points < min_num_overlapping_points)) {
    stop(paste0(
      "calc_extreme_shifts():\nafter applying stretch factor:", stretch, " to ", transformed.timecourse, ", none of the considered shifts have ",
      "min_num_overlapping_points (", min_num_overlapping_points, ") overlapping timepoints with the other timecourse!\n",
      "maybe try a smaller stretch, and double check you're applying it to the correct timecourse."
    ))
  }
  neg_extreme <- min(neg_extreme_candidate[num_overlapping_points >= min_num_overlapping_points])

  num_overlapping_points <- sapply(pos_extreme_candidates, FUN = calc_num_overlapping_points, original = original)
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
