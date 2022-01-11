.onLoad <- function(libname, pkgname) {
  options(cli.progress_show_after = 0)
  options(cli.progress_clear = FALSE)
}

#' Flag data to transform time points which overlap to reference data timecourse
#'
#' @noRd
get_compared_timepoints <- function(data,
                                    accession_data_to_transform,
                                    accession_data_ref) {
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

  # Use these to get all the Brassica times which used
  data$is_compared[(data$accession == accession_data_ref & (data$shifted_time >= min_data_ref & data$shifted_time <= max_data_ref))] <- TRUE

  return(data)
}

#' Calculate prediction of reference data expression value
#'
#' @noRd
interpolate_data_ref_comparison_expression <- function(data_to_transform_time, data_ref_dt) {
  # Suppress "no visible binding for global variable" note
  shifted_time <- NULL

  # Calculate diff
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

#' Validate names
#'
#' @noRd
match_names <- function(x, lookup) {
  unmatched <- x[-grep(paste(lookup, collapse = "$|"), x)]
  if (length(unmatched) > 0) {
    stop("Valid names are ", paste(lookup, collapse = ", "))
  }
}

#' Get approximate stretch factor
#'
#' @description
#' `get_approximate_stretch()` is a function to get a stretch factor estimation given input data. This function will take the time point ranges of both reference and query data and compare them to estimate the stretch factor.
#'
#' @param input_df Input data frame contains all replicates of gene expression in each genotype at each time point.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_ref Accession name of reference data.
#'
#' @return This function returns an estimation of a stretch factor for registering the data.
#' @export
get_approximate_stretch <- function(input_df, accession_data_to_transform, accession_data_ref) {
  # Suppress "no visible binding for global variable" note
  accession <- NULL

  # Calculate approximate stretch factor
  deltas <- input_df %>%
    dplyr::group_by(.data$accession) %>%
    dplyr::summarise(
      min_timepoint = min(.data$timepoint),
      max_timepoint = max(.data$timepoint),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      delta_timepoint = .data$max_timepoint - .data$min_timepoint
    )

  delta_ref <- deltas %>%
    dplyr::filter(accession == accession_data_ref) %>%
    dplyr::pull(.data$delta_timepoint)

  delta_to_transform <- deltas %>%
    dplyr::filter(accession == accession_data_to_transform) %>%
    dplyr::pull(.data$delta_timepoint)

  stretch_factor <- delta_ref / delta_to_transform

  return(stretch_factor)
}
