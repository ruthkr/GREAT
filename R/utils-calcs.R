#' Calculate Akaike's ‘An Information Criterion’
#'
#' @noRd
calc_AIC <- function(logL, num_params) {
  return((-2 * logL) + 2 * num_params)
}

#' Calculate Bayesian Information Criterion
#'
#' @noRd
calc_BIC <- function(logL, num_params, num_obs) {
  return((-2 * logL) + log(num_obs) * num_params)
}

#' Calculate mean score of observed expression data and expected value
#'
#' @noRd
calc_score <- function(data_to_transform_expression, data_ref_expression) {
  # if don't regularise / penalise for shift applied, then like uniform prior on it.
  # maybe should penalise for comparing fewer timepoints?
  # divide by the data_to_transform_expression as filtered already to make sure expressed
  score <- (data_to_transform_expression - data_ref_expression)**2 # / abs(data_to_transform_expression)

  return(mean(score))
}

#' Get the largest time point of reference data
#'
#' @noRd
max_is_compared_to_data_to_transform <- function(data_to_transform_time, data_ref) {
  # If using for rep data, then repeats of the same points screws it up
  data_ref <- unique(subset(data_ref, select = c("timepoint", "shifted_time")))

  # Return the reference data dt shifted timepoint which is greater than, or equal to the data to transform time.
  data_ref$diff <- data_ref$shifted_time - data_to_transform_time
  candidates <- data_ref[data_ref$diff >= 0, ]
  data_ref_max_time <- candidates$shifted_time[candidates$diff == min(candidates$diff)]

  return(data_ref_max_time)
}

#' Get the smallest time point of reference data
#'
#' @noRd
min_is_compared_to_data_to_transform <- function(data_to_transform_time, data_ref) {
  # If using for rep data, then repeats of the same points screws it up
  data_ref <- unique(subset(data_ref, select = c("timepoint", "shifted_time")))
  data_ref$diff <- data_ref$shifted_time - data_to_transform_time
  candidates <- data_ref[data_ref$diff <= 0, ]
  data_ref_min_time <- candidates$shifted_time[candidates$diff == max(candidates$diff)]

  return(data_ref_min_time)
}

#' Calculate the number of overlapping points
#'
#' @noRd
calc_num_overlapping_points <- function(shift, data, accession_data_to_transform, accession_data_ref) {
  # Suppress "no visible binding for global variable" note
  num.compared <- NULL
  is_compared <- NULL
  accession <- NULL

  # Calculations
  data$shifted_time[data$accession == accession_data_to_transform] <- data$delta_time[data$accession == accession_data_to_transform] + shift
  data <- get_compared_timepoints(
    data,
    accession_data_to_transform,
    accession_data_ref
  )
  data[, num.compared := sum(is_compared), by = .(accession)]

  return(min(data$num.compared))
}
