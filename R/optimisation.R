#' Optimise registration parameters with Simulated Annealing
#'
#' @param data Input data frame containing all replicates of gene expression for a single genotype at each time point.
#' @param overlapping_percent Number of minimum overlapping time points. Shifts will be only considered if it leaves at least these many overlapping points after applying the registration function.
#' @param optimisation_config List with arguments to modify the optimisation configuration.
#' @param optimise_fun Optimisation function to use. Can be \code{optimise_using_nm} or \code{optimise_using_nm}.
#'
#' @noRd
optimise <- function(data,
                     stretches = NA,
                     shifts = NA,
                     overlapping_percent = 0.5,
                     optimisation_config,
                     optimise_fun) {
  # Calculate boundary box and initial guess
  if (all(is.na(stretches), is.na(shifts))) {
    space_lims <- get_search_space_limits(data, overlapping_percent)
  } else {
    space_lims <- get_search_space_limits_from_params(stretches, shifts)
  }

  # Run optimisation
  optimised_params <- optimise_fun(
    data,
    optimisation_config,
    overlapping_percent,
    space_lims
  )

  return(optimised_params)
}

#' Objective loss function for Simulated Annealing
#'
#' @noRd
objective_fun <- function(data, stretch, shift, overlapping_percent, maximize = TRUE) {
  # Define objective function factor
  factor <- ifelse(maximize, 1, -1)

  tryCatch(
    {
      # Apply registration
      all_data_reg <- apply_registration(data, stretch, shift)

      # Check if overlapping condition is upheld
      if (calc_overlapping_percent(all_data_reg) < overlapping_percent) stop()

      # Calculate loglik
      loglik_combined <- factor * calc_loglik_H1(all_data_reg)
      return(loglik_combined)
    },
    error = function(error_message) {
      loglik_combined <- -factor * 999
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

#' Optimise stretch and shift using Simulated Annealing
#'
#' @noRd
optimise_using_sa <- function(data,
                              optimisation_config,
                              overlapping_percent,
                              space_lims) {
  # Parse initial and limit parameters
  stretch_init <- space_lims$stretch_init
  shift_init <- space_lims$shift_init
  stretch_lower <- space_lims$stretch_lower
  stretch_upper <- space_lims$stretch_upper
  shift_lower <- space_lims$shift_lower
  shift_upper <- space_lims$shift_upper

  # Optimisation parameters
  # TODO: Explore best default
  num_iterations <- optimisation_config$num_iterations
  num_inner_loop_iter <- optimisation_config$num_fun_evals

  # Calculate cooling schedule
  t0 <- 1000
  t_min <- 0.1
  r_cooling <- (t_min / t0)^(1 / num_iterations)

  # Perform SA using {optimization}
  optimised_params <- optimization::optim_sa(
    fun = function(x) objective_fun(data, x[1], x[2], overlapping_percent),
    maximization = TRUE,
    start = c(stretch_init, shift_init),
    trace = FALSE,
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

  # Results object
  params_list <- list(
    stretch = optimised_params$par[1],
    shift = optimised_params$par[2],
    loglik_score = optimised_params$function_value
  )

  return(params_list)
}


#' Optimise stretch and shift using Nelder-Mead
#'
#' @noRd
optimise_using_nm <- function(data,
                              optimisation_config,
                              overlapping_percent,
                              space_lims) {
  # Parse initial and limit parameters
  stretch_init <- space_lims$stretch_init
  shift_init <- space_lims$shift_init
  stretch_lower <- space_lims$stretch_lower
  stretch_upper <- space_lims$stretch_upper
  shift_lower <- space_lims$shift_lower
  shift_upper <- space_lims$shift_upper

  # Define data as it object required by optim::sa()
  fmsfundata <- structure(
    list(data = data),
    class = "optimbase.functionargs"
  )

  loglik_score_nm <- function(x = NULL, index = NULL, fmsfundata = NULL) {
    stretch <- x[1]
    shift <- x[2]

    f <- objective_fun(
      fmsfundata$data,
      stretch,
      shift,
      overlapping_percent,
      maximize = FALSE
    )

    varargout <- list(
      f = f,
      index = index,
      this = list(costfargument = fmsfundata)
    )

    return(varargout)
  }

  # Optimisation parameters
  # TODO: Explore best default
  num_iterations <- optimisation_config$num_iterations
  max_fun_evals <- optimisation_config$num_fun_evals

  # Start process optimisation
  x0 <- matrix(c(stretch_init, shift_init), ncol = 1)
  nm <- neldermead::neldermead()
  nm <- neldermead::neldermead.set(nm, "numberofvariables", 2)
  nm <- neldermead::neldermead.set(nm, "function", loglik_score_nm)
  nm <- neldermead::neldermead.set(nm, "x0", x0)
  nm <- neldermead::neldermead.set(nm, "costfargument", fmsfundata)
  nm <- neldermead::neldermead.set(nm, "maxiter", num_iterations)
  nm <- neldermead::neldermead.set(nm, "maxfunevals", max_fun_evals)
  nm <- neldermead::neldermead.set(nm, "method", "box")
  nm <- neldermead::neldermead.set(nm, "storehistory", FALSE)
  nm <- neldermead::neldermead.set(nm, "boundsmin", c(stretch_lower, shift_lower))
  nm <- neldermead::neldermead.set(nm, "boundsmax", c(stretch_upper, shift_upper))
  nm <- neldermead::neldermead.search(this = nm)

  # Parse simplex at optimal point
  simplex_obj <- unlist(nm$simplexopt)
  vertices <- nm$simplexopt$nbve
  simplex_vars <- grep("^x|^fv", names(simplex_obj), value = TRUE)

  optimised_params <- data.table::as.data.table(
    matrix(simplex_obj[simplex_vars], nrow = vertices)
  )[vertices, ]

  # Results object
  params_list <- list(
    stretch = optimised_params$V1,
    shift = optimised_params$V2,
    loglik_score = -optimised_params$V3
  )

  return(params_list)
}

#' Optimise stretch and shift using L-BFGS-B
#'
#' @noRd
optimise_using_lbfgsb <- function(data,
                                  optimisation_config = NULL,
                                  overlapping_percent,
                                  space_lims) {

  # Parse initial and limit parameters
  stretch_init <- space_lims$stretch_init
  shift_init <- space_lims$shift_init
  stretch_lower <- space_lims$stretch_lower
  stretch_upper <- space_lims$stretch_upper
  shift_lower <- space_lims$shift_lower
  shift_upper <- space_lims$shift_upper

  loglik_score_lbfgsb <- function(theta, data_input, overlapping_percent) {
    stretch <- theta[1]
    shift <- theta[2]

    loglik_score <- objective_fun(
      data_input,
      stretch,
      shift,
      overlapping_percent,
      maximize = FALSE
    )

    return(loglik_score)
  }

  results <- stats::optim(
    par = c(stretch_init, shift_init),
    fn = loglik_score_lbfgsb,
    data_input = data,
    overlapping_percent = overlapping_percent,
    method = "L-BFGS-B",
    lower = c(stretch_lower, shift_lower),
    upper = c(stretch_upper, shift_upper)
  )

  # Results object
  params_list <- list(
    stretch = results$par[1],
    shift = results$par[2],
    loglik_score = -results$value
  )

  return(params_list)
}


