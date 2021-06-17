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

# ara.expression <- ara.compared$mean.cpm
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
  bra.dt <- unique(subset(bra.dt, select = c("timepoint", "shifted.time")))

  # return the bra dt shifted timepoint which is greater than, or equal to the Ara time.
  bra.dt$diff <- bra.dt$shifted.time - arabidopsis.time
  candidates <- bra.dt[bra.dt$diff >= 0, ]
  bra.max.time <- candidates$shifted.time[candidates$diff == min(candidates$diff)]

  return(bra.max.time)

  # data.table::setorder(bra.dt, diff)
  # nearest.points <- bra.dt[1:2,]
  # data.table::setorder(nearest.points, shifted.time)
  # return(nearest.points$shifted.time[2])
}

# arabidopsis.time=ara.min
#' @export
min_is_compared_to_arabidopsis <- function(arabidopsis.time, bra.dt) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # return the smallest brassica time which is used in comparison t the arabidopsis time
  # the biggest one smaller than or equal to the arabidopsis time

  # if using for rep data, then repeats of the same points screws it up
  bra.dt <- unique(subset(bra.dt, select = c("timepoint", "shifted.time")))

  bra.dt$diff <- bra.dt$shifted.time - arabidopsis.time
  candidates <- bra.dt[bra.dt$diff <= 0, ]
  bra.min.time <- candidates$shifted.time[candidates$diff == max(candidates$diff)]

  return(bra.min.time)
}


#' @export
calc_num_overlapping_points <- function(shift, original) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # calculate the number of overlapping points for the species with the fewer overlapping points if the current "shift" is
  # applied to the col0 delta timepoints.
  original$shifted.time[original$accession == "Col0"] <- original$delta_time[original$accession == "Col0"] + shift
  original <- get_compared_timepoints(original)
  original[, num.compared := sum(is.compared), by = .(accession)]

  return(min(original$num.compared))
}
