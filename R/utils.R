#' @export
get_compared_timepoints <- function(test) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # flag the arabidopsis timepoints which overlap the brassica timecourse, and so will be compared
  bra.min <- min(test$shifted_time[test$accession=='Ro18'])
  bra.max <- max(test$shifted_time[test$accession=='Ro18'])

  # get the arabidopsis times which used
  test$is.compared <- FALSE
  test$is.compared[(test$accession=='Col0' & (test$shifted_time >= bra.min & test$shifted_time <=bra.max))] <- TRUE

  # get the extreme brassica times which used - bigger or equal than Ara max, and smaller or equal than Ara min, because have to project
  #  Ara onto Bra
  ara.max <- max(test$shifted_time[test$accession=='Col0' & test$is.compared==TRUE])
  ara.min <- min(test$shifted_time[test$accession=='Col0' & test$is.compared==TRUE])
  bra.max <- max_is_compared_to_arabidopsis(ara.max, test[test$accession=='Ro18', ])
  bra.min <- min_is_compared_to_arabidopsis(ara.min, test[test$accession=='Ro18', ])

  # use these to get all the brassica times which used
  test$is.compared[(test$accession=='Ro18' & (test$shifted_time >= bra.min & test$shifted_time <=bra.max))] <- TRUE

  return(test)
}


# arabidopsis.time <- 2
# bra.dt <- ara.df
#' @export
interpolate_brassica_comparison_expression <- function(arabidopsis.time, bra.dt) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # arabidopsis time is outside of the range of the bra.dt shifted timepoints
  bra.dt$diff <- bra.dt$shifted_time - arabidopsis.time

  # if outside of comparible range (time is smaller than all bra.dt time or bigger than all)
  if (all(bra.dt$diff > 0) | all(bra.dt$diff < 0) ) {
    return(NA)
  }

  # otherwise,  cut down brassica observations to the two nearest timepoints to the arabidopsis time
  bra.dt$diff <- abs(bra.dt$shifted_time - arabidopsis.time)
  data.table::setorder(bra.dt, diff)
  nearest.points <- bra.dt[1:2,]

  # linearly interpolate between these points to estimate the comparison expression value
  data.table::setorder(nearest.points, shifted_time) # so [1] is earlier time
  time.diff <- nearest.points$shifted_time[2] - nearest.points$shifted_time[1] #
  expression.diff <- nearest.points$mean_cpm[2] - nearest.points$mean_cpm[1]
  grad <- expression.diff / time.diff
  pred.expression <- nearest.points$mean_cpm[1] + (nearest.points$diff[1]) * grad

  return(pred.expression)
}
