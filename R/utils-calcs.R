#' Scaling value
#'
#' `my_scale` is a function to get a scaled value of an input by dividing the maximum value.
#'
#' @param input Input data.
#'
#' @return Scaled value.
#' @export
my_scale <- function(input) {
  return(input / max(input))
}

#' Calculate Akaike's ‘An Information Criterion’
#'
#' `calc_AIC` is a generic function calculating Akaike's ‘An Information Criterion’ for one or several fitted model objects for which a log-likelihood value can be obtained.
#'
#' @param logL Loglikelihoods value obtained by stats::logLik().
#' @param num_params Number of parameters in the fitted model.
#'
#' @return AIC criterion value.
#' @export
calc_AIC <- function(logL, num_params) {
  return((-2 * logL) + 2 * num_params)
}

#' Calculate Bayesian Information Criterion
#'
#' `calc_BIC` is a function calculating Bayesian Information Criterion, which is a special case of AIC when k = log(n), where n being the number of observations.
#'
#' @param logL Loglikelihoods value obtained by stats::logLik().
#' @param num_params Number of parameters in the fitted model.
#' @param num_obs Number of observations.
#'
#' @return AIC criterion value.
#' @export
calc_BIC <- function(logL, num_params, num_obs) {
  return((-2 * logL) + log(num_obs) * num_params)
}

#' Calculate mean score of observed expression data and expected value
#'
#' `calc_score` is a function to calculate a mean of a score of difference of observed expression data and expected expression value, by executing score = (sum(observed-expected)**2). The mean is taken to compare the variable number of datapoints.
#'
#' @param data_to_transform_expression Input expression of data_to_transform.
#' @param data_ref_expression Input expression of data_ref.
#'
#' @return Mean of score value.
#' @export
calc_score <- function(data_to_transform_expression, data_ref_expression) {
  # if don't regularise / penalise for shift applied, then like uniform prior on it.
  # maybe should penalise for comparing fewer timepoints?
  # divide by the data_to_transform_expression as filtered already to make sure expressed
  score <- (data_to_transform_expression - data_ref_expression)**2 # / abs(data_to_transform_expression)

  return(mean(score))
}

#' Get the largest time point of reference data
#'
#' `max_is_compared_to_data_to_transform` is used to get the largest time reference data which is used in comparison t the data_to_transform time.
#'
#' @param data_to_transform_time Maximum time points of candidate data to transform.
#' @param data_ref A data frame containing reference data.
#'
#' @export
max_is_compared_to_data_to_transform <- function(data_to_transform_time,
                                                 data_ref) {
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
#' `min_is_compared_to_data_to_transform` is used to get the smallest time reference data which is used in comparison t the data_to_transform time.
#'
#' @param data_to_transform_time Minimum time points of candidate data to transform.
#' @param data_ref A data frame containing reference data.
#'
#' @export
min_is_compared_to_data_to_transform <- function(data_to_transform_time,
                                                 data_ref) {
  # If using for rep data, then repeats of the same points screws it up
  data_ref <- unique(subset(data_ref, select = c("timepoint", "shifted_time")))
  data_ref$diff <- data_ref$shifted_time - data_to_transform_time
  candidates <- data_ref[data_ref$diff <= 0, ]
  data_ref_min_time <- candidates$shifted_time[candidates$diff == max(candidates$diff)]

  return(data_ref_min_time)
}

#' Calculate the number of overlapping points
#'
#' `calc_num_overlapping_points` is used to calculate the number of overlapping points for the species with the fewer overlapping points if the current "shift" is applied to the data to alined delta timepoints.
#'
#' @param data Input data.
#' @param shift Current shift value
#' @param accession_data_to_transform Accession name of data which will be transformed.
#'
#' @export
calc_num_overlapping_points <- function(shift, data, accession_data_to_transform = "Col0") {
  data$shifted_time[data$accession == accession_data_to_transform] <- data$delta_time[data$accession == accession_data_to_transform] + shift
  data <- get_compared_timepoints(data)
  data[, num.compared := sum(is_compared), by = .(accession)]

  return(min(data$num.compared))
}
