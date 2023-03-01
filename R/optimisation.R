#' Optimise registration parameters with Simulated Annealing
#'
#' @param data Input data frame containing all replicates of gene expression for a single genotype at each time point.
#' @param overlapping_percent Number of minimum overlapping time points. Shifts will be only considered if it leaves at least these many overlapping points after applying the registration function.
#' @param optimisation_config List with optional arguments to modify the Simulated Annealing optimisation.
#'
#' @noRd
optimise <- function(data,
                     stretches = NA,
                     shifts = NA,
                     overlapping_percent = 0.5,
                     optimisation_config) {
  # Calculate boundary box and initial guess
  if (all(is.na(stretches), is.na(shifts))) {
    space_lims <- get_search_space_limits(data, overlapping_percent)
  } else {
    space_lims <- get_search_space_limits_from_params(stretches, shifts)
  }

  # Parse initial and limit parameters
  stretch_init <- space_lims$stretch_init
  shift_init <- space_lims$shift_init
  stretch_lower <- space_lims$stretch_lower
  stretch_upper <- space_lims$stretch_upper
  shift_lower <- space_lims$shift_lower
  shift_upper <- space_lims$shift_upper

  # Calculate cooling schedule
  t0 <- 1000
  t_min <- 0.1
  r_cooling <- (t_min / t0)^(1 / optimisation_config$num_iterations)
  # TODO: Explore best default
  num_inner_loop_iter <- 100

  # Perform SA using {optimization}
  optimised_params <- optimization::optim_sa(
    fun = function(x) objective_fun(data, x[1], x[2], overlapping_percent),
    maximization = TRUE,
    start = c(stretch_init, shift_init),
    trace = TRUE,
    lower = c(stretch_lower, shift_lower),
    upper = c(stretch_upper, shift_upper),
    control = list(
      t0 = t0,
      t_min = t_min,
      nlimit = num_inner_loop_iter,
      r = r_cooling,
      rf = 3,
      dyn_rf = FALSE
    )
  )

  return(optimised_params)
}

#' Objective loss function for Simulated Annealing
#'
#' @noRd
objective_fun <- function(data, stretch, shift, overlapping_percent) {
  tryCatch(
    {
      # Apply registration
      all_data_reg <- apply_registration(data, stretch, shift)

      # Check if overlapping condition is upheld
      if (calc_overlapping_percent(all_data_reg) < overlapping_percent) stop()

      # Calculate loglik
      loglik_combined <- calc_loglik_H1(all_data_reg)
      return(loglik_combined)
    },
    error = function(error_message) {
      loglik_combined <- -999
      return(loglik_combined)
    }
  )
}

#' Calculate limits of the search space for Simulated Annealing
#'
#' @noRd
get_search_space_limits <- function(data, overlapping_percent = 0.5) {
  # Suppress "no visible binding for global variable" note
  accession <- NULL
  timepoint <- NULL

  # Initial stretch limits
  stretch_init <- get_approximate_stretch(data)
  stretch_lower <- 0.5 * stretch_init
  stretch_upper <- 1.5 * stretch_init

  # Extract time point ranges
  timepoints_ref <- unique(data[accession == "ref", timepoint])
  timepoints_query <- unique(data[accession == "query", timepoint])

  # Calculate time point ranges
  range_ref <- diff(range(timepoints_ref))
  range_query <- diff(range(timepoints_query))
  range_query_max_stretch <- stretch_upper * range_query

  # Calculate minimum and maximum timepoints in which the curves overlap
  min_timepoint <- min(timepoints_ref) + overlapping_percent * range_ref - range_query_max_stretch
  max_timepoint <- max(timepoints_ref) - overlapping_percent * range_ref + range_query_max_stretch

  # Calculate shift limits
  shift_lower <- min_timepoint - min(timepoints_query)
  shift_upper <- (max_timepoint - range_query_max_stretch) - min(timepoints_query)

  # Calculate initial shift value (zero if possible)
  shift_init <- 0
  if (shift_init < shift_lower | shift_init > shift_upper) {
    shift_init <- mean(c(shift_lower, shift_upper))
  }

  # Results object
  results_list <- list(
    stretch_init = stretch_init,
    stretch_lower = stretch_lower,
    stretch_upper = stretch_upper,
    shift_init = shift_init,
    shift_lower = shift_lower,
    shift_upper = shift_upper
  )

  return(results_list)
}

#' Calculate limits of the search space for Simulated Annealing from provided registration parameters
#'
#' @noRd
get_search_space_limits_from_params <- function(stretches, shifts) {
  # Initial stretch limits
  stretch_lower <- min(stretches)
  stretch_upper <- max(stretches)
  stretch_init <- mean(stretches)

  # Calculate shift limits
  shift_lower <- min(shifts)
  shift_upper <- max(shifts)

  # Calculate initial shift value (zero if possible)
  shift_init <- 0
  if (shift_init < shift_lower | shift_init > shift_upper) {
    shift_init <- mean(c(shift_lower, shift_upper))
  }

  # Results object
  results_list <- list(
    stretch_init = stretch_init,
    stretch_lower = stretch_lower,
    stretch_upper = stretch_upper,
    shift_init = shift_init,
    shift_lower = shift_lower,
    shift_upper = shift_upper
  )

  return(results_list)
}

#' Calculate overlapping percentage between reference and query data time point ranges
#'
#' @noRd
calc_overlapping_percent <- function(data) {
  # Suppress "no visible binding for global variable" note
  accession <- NULL
  timepoint <- NULL

  # Extract time point ranges
  range_ref <- range(unique(data[accession == "ref", timepoint]))
  range_query <- range(unique(data[accession == "query", timepoint]))

  if (all(range_ref[2] >= range_query[2], range_ref[1] <= range_query[1])) {
    # Query is fully contained on reference
    overlapping_percent <- 1
  } else {
    # Calculate overlapping percent over reference
    overlap <- min(c(range_ref[2], range_query[2])) - max(c(range_ref[1], range_query[1]))
    overlapping_percent <- overlap / diff(range_ref)
  }

  return(overlapping_percent)
}
