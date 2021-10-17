.onLoad <- function(libname, pkgname) {
  options(cli.progress_show_after = 0)
  options(cli.progress_clear = FALSE)
}

#' Flag data to transform time points which overlap to reference data timecourse
#'
#' `get_compared_timepoints` flags data to transform time points which overlap to reference data timecourse by comparing each data to transform time points to minimum and maximum value of reference data.
#'
#' @param data Input data containing both data to transform and reference data.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_ref Accession name of reference data.
get_compared_timepoints <- function(data,
                                    accession_data_to_transform = "Col0",
                                    accession_data_ref = "Ro18") {
  # Filter reference data from the whole dataset
  data_ref <- data$shifted_time[data$accession == accession_data_ref]

  min_data_ref <- min(data_ref)
  max_data_ref <- max(data_ref)

  # Get time points of data to transform which are used
  data$is_compared <- FALSE
  data$is_compared[(data$accession == accession_data_to_transform & (data$shifted_time >= min_data_ref & data$shifted_time <= max_data_ref))] <- TRUE

  # Get the extreme reference data times which used - bigger or equal than max of data to transform, and smaller or equal than  min data to transform, because have to project data to transform onto reference data
  max_data_to_transform <- max(data$shifted_time[data$accession == accession_data_to_transform & data$is_compared == TRUE])
  min_data_to_transform <- min(data$shifted_time[data$accession == accession_data_to_transform & data$is_compared == TRUE])
  max_data_ref <- max_is_compared_to_data_to_transform(
    data_to_transform_time = max_data_to_transform,
    data_ref = data[data$accession == accession_data_ref, ]
  )
  min_data_ref <- min_is_compared_to_data_to_transform(
    data_to_transform_time = min_data_to_transform,
    data_ref = data[data$accession == accession_data_ref, ]
  )

  # use these to get all the brassica times which used
  data$is_compared[(data$accession == accession_data_ref & (data$shifted_time >= min_data_ref & data$shifted_time <= max_data_ref))] <- TRUE

  return(data)
}

#' Calculate prediction of reference data expression value
#'
#' @param data_to_transform_time Input time from data to transform.
#' @param data_ref_dt Input reference data data frame.
#'
#' @return Expression prediction.
interpolate_data_ref_comparison_expression <- function(data_to_transform_time,
                                                       data_ref_dt) {
  data_ref_dt$diff <- data_ref_dt$shifted_time - data_to_transform_time

  # If outside of comparable range (time is smaller than all data_ref_dt time or bigger than all)
  if (all(data_ref_dt$diff > 0) | all(data_ref_dt$diff < 0)) {
    return(NA)
  }

  # Otherwise, cut down data_ref observations to the two nearest timepoints to the data_to_transform time
  data_ref_dt$diff <- abs(data_ref_dt$shifted_time - data_to_transform_time)
  data.table::setorder(data_ref_dt, diff)
  nearest.points <- data_ref_dt[1:2, ]

  # Linearly interpolate between these points to estimate the comparison expression value
  data.table::setorder(nearest.points, shifted_time) # so [1] is earlier time
  time.diff <- nearest.points$shifted_time[2] - nearest.points$shifted_time[1]
  expression.diff <- nearest.points$expression_value[2] - nearest.points$expression_value[1]
  grad <- expression.diff / time.diff
  pred.expression <- nearest.points$expression_value[1] + (nearest.points$diff[1]) * grad

  return(pred.expression)
}
