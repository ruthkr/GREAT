my.scale <- function(v) {
  return(v / max(v))
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

# arabidopsis.time <- ara.max
# arabidopsis.time <- ara.max
# bra.dt <- test[test$accession=='Ro18',]
#' @export
max_is_compared_to_arabidopsis <- function(arabidopsis.time, bra.dt) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # return the largest brassica time which is used in comparison t the arabidopsis time
  # the smallest one greater to or equal to arabidopsis time

  # if using for rep data, then repeats of the same points screws it up
  bra.dt <- unique(subset(bra.dt, select = c("timepoint", "shifted_time")))

  # return the bra dt shifted timepoint which is greater than, or equal to the Ara time.
  bra.dt$diff <- bra.dt$shifted_time - arabidopsis.time
  candidates <- bra.dt[bra.dt$diff >= 0, ]
  bra.max.time <- candidates$shifted_time[candidates$diff == min(candidates$diff)]

  return(bra.max.time)

  # data.table::setorder(bra.dt, diff)
  # nearest.points <- bra.dt[1:2,]
  # data.table::setorder(nearest.points, shifted_time)
  # return(nearest.points$shifted_time[2])
}

# arabidopsis.time=ara.min
#' @export
min_is_compared_to_arabidopsis <- function(arabidopsis.time, bra.dt) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # return the smallest brassica time which is used in comparison t the arabidopsis time
  # the biggest one smaller than or equal to the arabidopsis time

  # if using for rep data, then repeats of the same points screws it up
  bra.dt <- unique(subset(bra.dt, select = c("timepoint", "shifted_time")))

  bra.dt$diff <- bra.dt$shifted_time - arabidopsis.time
  candidates <- bra.dt[bra.dt$diff <= 0, ]
  bra.min.time <- candidates$shifted_time[candidates$diff == max(candidates$diff)]

  return(bra.min.time)
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
