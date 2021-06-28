#' Flag data to align time points which overlap to data target timecourse
#'
#' `get_compared_timepoints` flags data to align time points which overlap to data target timecourse by comparing each data to align time points to minimum and maximum value of data target.
#'
#' @param data Input data containing both data to align and data target.
#' @param accession_data_to_align Accession name of data which will be aligned.
#' @param accession_data_target Accession name of data target.
#'
#' @export
get_compared_timepoints <- function(data,
                                    accession_data_to_align = "Col0",
                                    accession_data_target = "Ro18") {

  # Filter data target from the whole dataset
  data_target <- data$shifted_time[data$accession == accession_data_target]

  min_data_target <- min(data_target)
  max_data_target <- max(data_target)

  # Get time points of data to align which are used
  data$is_compared <- FALSE
  data$is_compared[(data$accession == accession_data_to_align & (data$shifted_time >= min_data_target & data$shifted_time <= max_data_target))] <- TRUE

  # Get the extreme data target times which used - bigger or equal than max of data to align, and smaller or equal than  min data to align, because have to project data to align onto data target
  max_data_to_align <- max(data$shifted_time[data$accession == accession_data_to_align & data$is_compared==TRUE])
  min_data_to_align <- min(data$shifted_time[data$accession == accession_data_to_align & data$is_compared==TRUE])
  max_data_target <- max_is_compared_to_data_to_align(max_data_to_align, data[data$accession == accession_data_target, ])
  min_data_target <- min_is_compared_to_data_to_align(min_data_to_align, data[data$accession == accession_data_target, ])

  # use these to get all the brassica times which used
  data$is_compared[(data$accession == accession_data_target & (data$shifted_time >= min_data_target & data$shifted_time <=max_data_target))] <- TRUE

  return(data)
}



#' Calculate prediction of data target expression value
#'
#' @param data_to_align_time Input time from data to align.
#' @param data_target_dt Input data target data frame.
#'
#' @return Expression prediction.
interpolate_data_target_comparison_expression <- function(data_to_align_time,
                                                          data_target_dt) {

  data_target_dt$diff <- data_target_dt$shifted_time - data_to_align_time

  # If outside of comparable range (time is smaller than all data_target_dt time or bigger than all)
  if (all(data_target_dt$diff > 0) | all(data_target_dt$diff < 0)) {
    return(NA)
  }

  # Otherwise, cut down data_target observations to the two nearest timepoints to the data_to_align time
  data_target_dt$diff <- abs(data_target_dt$shifted_time - data_to_align_time)
  data.table::setorder(data_target_dt, diff)
  nearest.points <- data_target_dt[1:2, ]

  # Linearly interpolate between these points to estimate the comparison expression value
  data.table::setorder(nearest.points, shifted_time) # so [1] is earlier time
  time.diff <- nearest.points$shifted_time[2] - nearest.points$shifted_time[1]
  expression.diff <- nearest.points$mean_cpm[2] - nearest.points$mean_cpm[1]
  grad <- expression.diff / time.diff
  pred.expression <- nearest.points$mean_cpm[1] + (nearest.points$diff[1]) * grad

  return(pred.expression)
}
