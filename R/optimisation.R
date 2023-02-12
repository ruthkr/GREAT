#' Optimise registration parameters with Simulated Annealing
#'
#' @param data Input data frame containing all replicates of gene expression for a single genotype at each time point.
#' @param overlapping_percent Number of minimum overlapping time points. Shifts will be only considered if it leaves at least these many overlapping points after applying the registration function.
#' @param optimisation_config List with optional arguments to modify the Simulated Annealing optimisation.
#'
#' @noRd
optimise <- function(data,
                     overlapping_percent = 0.5,
                     optimisation_config) {
  # Calculate boundary box and initial guess
  space_lims <- get_search_space_limits(data, overlapping_percent)

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
    fun = function(x) objective_fun(data, x[1], x[2]),
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
objective_fun <- function(data, stretch, shift) {
  tryCatch(
    {
      all_data_reg <- apply_registration(data, stretch, shift)
      BIC_combined <- compare_dynamics_H1(all_data_reg)
      return(BIC_combined)
    },
    error = function(error_message) {
      BIC_combined <- 999
      return(BIC_combined)
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
  stretch_lower <- 0.75 * stretch_init
  stretch_upper <- 1.25 * stretch_init

  # Extract timepoint ranges
  timepoints_ref <- unique(data[accession == "ref", timepoint])
  timepoints_query <- unique(data[accession == "query", timepoint])

  # Calculate timepoint ranges
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
