#' Flag data to transform time points which overlap to data fix timecourse
#'
#' `get_compared_timepoints` flags data to transform time points which overlap to data fix timecourse by comparing each data to transform time points to minimum and maximum value of data fix.
#'
#' @param data Input data containing both data to transform and data fix.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_fix Accession name of data fix.
#'
#' @export
get_compared_timepoints <- function(data,
                                    accession_data_to_transform = "Col0",
                                    accession_data_fix = "Ro18") {
  # Filter data fix from the whole dataset
  data_fix <- data$shifted_time[data$accession == accession_data_fix]

  min_data_fix <- min(data_fix)
  max_data_fix <- max(data_fix)

  # Get time points of data to transform which are used
  data$is_compared <- FALSE
  data$is_compared[(data$accession == accession_data_to_transform & (data$shifted_time >= min_data_fix & data$shifted_time <= max_data_fix))] <- TRUE

  # Get the extreme data fix times which used - bigger or equal than max of data to transform, and smaller or equal than  min data to transform, because have to project data to transform onto data fix
  max_data_to_transform <- max(data$shifted_time[data$accession == accession_data_to_transform & data$is_compared == TRUE])
  min_data_to_transform <- min(data$shifted_time[data$accession == accession_data_to_transform & data$is_compared == TRUE])
  max_data_fix <- max_is_compared_to_data_to_transform(
    data_to_transform_time = max_data_to_transform,
    data_fix = data[data$accession == accession_data_fix, ]
  )
  min_data_fix <- min_is_compared_to_data_to_transform(
    data_to_transform_time = min_data_to_transform,
    data_fix = data[data$accession == accession_data_fix, ]
  )

  # use these to get all the brassica times which used
  data$is_compared[(data$accession == accession_data_fix & (data$shifted_time >= min_data_fix & data$shifted_time <= max_data_fix))] <- TRUE

  return(data)
}

#' Calculate prediction of data fix expression value
#'
#' @param data_to_transform_time Input time from data to transform.
#' @param data_fix_dt Input data fix data frame.
#'
#' @return Expression prediction.
interpolate_data_fix_comparison_expression <- function(data_to_transform_time,
                                                       data_fix_dt) {
  data_fix_dt$diff <- data_fix_dt$shifted_time - data_to_transform_time

  # If outside of comparable range (time is smaller than all data_fix_dt time or bigger than all)
  if (all(data_fix_dt$diff > 0) | all(data_fix_dt$diff < 0)) {
    return(NA)
  }

  # Otherwise, cut down data_fix observations to the two nearest timepoints to the data_to_transform time
  data_fix_dt$diff <- abs(data_fix_dt$shifted_time - data_to_transform_time)
  data.table::setorder(data_fix_dt, diff)
  nearest.points <- data_fix_dt[1:2, ]

  # Linearly interpolate between these points to estimate the comparison expression value
  data.table::setorder(nearest.points, shifted_time) # so [1] is earlier time
  time.diff <- nearest.points$shifted_time[2] - nearest.points$shifted_time[1]
  expression.diff <- nearest.points$mean_cpm[2] - nearest.points$mean_cpm[1]
  grad <- expression.diff / time.diff
  pred.expression <- nearest.points$mean_cpm[1] + (nearest.points$diff[1]) * grad

  return(pred.expression)
}
