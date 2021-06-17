# test <- mean_df
# stretch_factor
# min_num_overlapping_points
# shift_extreme
#' @export
get_extreme_shifts_for_all <- function(test, stretch_factor, min_num_overlapping_points, shift_extreme) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # wrapper for calc_extreme_shifts to be able to move it out of the loop so don't calculate for every gene.

  #min_num_overlapping_points <- 5 # bound the extreme allowed shifts, such than at least this many timepoints are being compared for both accessions.
  # cut data.table to a single gene
  curr_sym <- unique(test$locus_name)[1]
  test <- test[test$locus_name==curr_sym, ]

  # transform timepoint to be time from first timepoint
  test[, delta_time:=timepoint - min(timepoint), by=.(accession)]
  # apply stretch_factor to the arabidopsis, leave the rapa as is
  test$delta_time[test$accession=='Col0'] <- test$delta_time[test$accession=='Col0']*stretch_factor

  # calculate min shift and max time shift, which still allows overlap of at least 5 times to be compared from whichever accession will be considering fewer timepoints from.
  # Shift is applied to the arabidopsis - so the 5th largest arabidopsis time is the biggest -ve shift can be applied
  # and the biggest shift which can be applied is to make the 5th smallest arabidopsis time == largest brassica time
  #data.table::setorder(test, delta_time)
  #min_shift <- min(test$delta_time[test$accession=='Ro18']) - test$delta_time[test$accession=='Col0'][length(test$delta_time[test$accession=='Col0'])-4]
  #max_shift <- max(test$delta_time[test$accession=='Ro18']) - test$delta_time[test$accession=='Col0'][5]

  M <- calc_extreme_shifts(test, min_num_overlapping_points, shift_extreme)
  return(M)
}

#' @export
calc_extreme_shifts <- function(test, min_num_overlapping_points, shift_extreme) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # calculate the minimum and maximum shifts can apply to Col-0 after the stretch transformation, whilst
  # preserving the criteria that at least min_num_overlapping_points are being compared from both accessions.

  original <- data.table::copy(test)
  original$shifted_time <- original$delta_time

  # print('line 1803')
  # print(original)

  # -ve extreme shift will be -1*exactly the difference between 1 of the stretched Col0 timepoints, and the smallest Ro18 timepoint
  # +ve extreme will be the difference between 1 of the col0 timepoints, and the maximum Ro18 timepoint
  neg_extreme_candidate <- -1*(original$delta_time[original$accession=='Col0'] - min(original$delta_time[original$accession=='Ro18']))
  pos_extreme_candidates <- max(original$delta_time[original$accession=='Ro18']) - original$delta_time[original$accession=='Col0']

  # of these candidates, find the most extreme values which mainting the required number of overlapping timepoints to be considered.
  num_overlapping_points <- sapply(neg_extreme_candidate, FUN=calc_num_overlapping_points, original=original)
  if (all(num_overlapping_points < min_num_overlapping_points)) {
    stop(paste0('calc_extreme_shifts():\nafter applying stretch factor:', stretch, ' to ', transformed.timecourse, ', none of the considered shifts have ',
                'min_num_overlapping_points (', min_num_overlapping_points, ') overlapping timepoints with the other timecourse!\n',
                "maybe try a smaller stretch, and double check you're applying it to the correct timecourse." ))
  }

  neg_extreme <- min(neg_extreme_candidate[num_overlapping_points >= min_num_overlapping_points])

  num_overlapping_points <- sapply(pos_extreme_candidates, FUN=calc_num_overlapping_points, original=original)
  pos_extreme <- max(pos_extreme_candidates[num_overlapping_points >= min_num_overlapping_points])

  # hard code maximum and minimum allowed shifts, as noticed spurious registrations when too extreme shifts
  # allowed
  if (neg_extreme < (-1*shift_extreme)) {
    neg_extreme <- -1 * shift_extreme
  }
  if (pos_extreme > 1*shift_extreme) {
    pos_extreme <- shift_extreme
  }

  return(list(neg_extreme, pos_extreme))
}
