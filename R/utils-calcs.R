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

# calculate the comparison stats
#' @export
calc.AIC <- function(logL, num.params) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  return((-2 * logL) + 2 * num.params)
}

#' @export
calc.BIC <- function(logL, num.params, num.obs) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  return((-2 * logL) + log(num.obs) * num.params)
}

# ara.expression <- ara.compared$mean_cpm
# bra.expression <- ara.compared$pred.bra.expression
#' @export
calc_score <- function(ara.expression, bra.expression) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # (sum(observed-expected)**2)
  # take mean, because going to be comparing variable number of datapoints.
  # if don't regularise / penalise for shift applied, then like uniform prior on it.
  # maybe should penalise for comparing fewer timepoints?
  # divide by the ara.expression as filtered already to make sure expressed

  d <- (ara.expression - bra.expression)**2 # / abs(ara.expression)
  return(mean(d))
}

#' Get the largest time point of data target
#'
#' `max_is_compared_to_data_to_align` is used to get the largest time data target which is used in comparison t the data_to_align time.
#'
#' @param data_to_align_time Maximum time points of candidate data to align.
#' @param data_target A data frame containing data target.
#'
#' @export
max_is_compared_to_data_to_align <- function(data_to_align_time,
                                             data_target) {

  # If using for rep data, then repeats of the same points screws it up
  data_target <- unique(subset(data_target, select = c("timepoint", "shifted_time")))

  # Return the data target dt shifted timepoint which is greater than, or equal to the data to align time.
  data_target$diff <- data_target$shifted_time - data_to_align_time
  candidates <- data_target[data_target$diff >= 0, ]
  data_target_max_time <- candidates$shifted_time[candidates$diff == min(candidates$diff)]

  return(data_target_max_time)

}


#' Get the smallest time point of data target
#'
#' `min_is_compared_to_data_to_align` is used to get the smallest time data target which is used in comparison t the data_to_align time.
#'
#' @param data_to_align_time Minimum time points of candidate data to align.
#' @param data_target A data frame containing data target.
#'
#' @export
min_is_compared_to_data_to_align <- function(data_to_align_time,
                                             data_target) {

  # If using for rep data, then repeats of the same points screws it up
  data_target <- unique(subset(data_target, select = c("timepoint", "shifted_time")))

  data_target$diff <- data_target$shifted_time - data_to_align_time
  candidates <- data_target[data_target$diff <= 0, ]
  data_target_min_time <- candidates$shifted_time[candidates$diff == max(candidates$diff)]

  return(data_target_min_time)
}


#' Calculate the number of overlapping points
#'
#' `calc_num_overlapping_points` is used to calculate the number of overlapping points for the species with the fewer overlapping points if the current "shift" is applied to the data to alined delta timepoints.
#'
#' @param data Input data.
#' @param shift Current shift value
#' @param accession_data_to_align Accession name of data which will be aligned.
#'
#' @export
calc_num_overlapping_points <- function(shift, data, accession_data_to_align = "Col0") {

  data$shifted_time[data$accession == accession_data_to_align] <- data$delta_time[data$accession == accession_data_to_align] + shift
  data <- get_compared_timepoints(data)

  data[, num.compared := sum(is_compared), by = .(accession)]

  return(min(data$num.compared))
}
