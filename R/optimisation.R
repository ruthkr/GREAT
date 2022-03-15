#' Calculate boundary box for Simulated Annealing
#'
#' @noRd
get_boundary_box <- function(input_df, accession_data_to_transform, accession_data_ref,
                             min_num_overlapping_points, expression_value_threshold, boundary_coverage) {
  # Stretch limits
  stretch_init <- get_approximate_stretch(
    input_df = input_df,
    accession_data_to_transform = accession_data_to_transform,
    accession_data_ref = accession_data_ref
  )

  stretch_lower <- 1
  stretch_upper <- 2 * ceiling(stretch_init) - stretch_lower

  # Shift limits
  mean_df <- get_mean_data(
    exp = input_df,
    expression_value_threshold = expression_value_threshold,
    accession_data_to_transform = accession_data_to_transform,
    is_data_normalised = FALSE
  )

  shift_upper <- get_extreme_shifts_for_all(
    mean_df,
    stretch_factor = 1.5,
    min_num_overlapping_points = min_num_overlapping_points,
    shift_extreme = 1000,
    accession_data_to_transform = accession_data_to_transform,
    accession_data_ref = accession_data_ref
  ) %>%
    .[[2]]

  shift_lower <- get_extreme_shifts_for_all(
    mean_df,
    stretch_factor = stretch_upper,
    min_num_overlapping_points = 4,
    shift_extreme = 1000,
    accession_data_to_transform = accession_data_to_transform,
    accession_data_ref = accession_data_ref
  )  %>%
    .[[1]]

  # Results object
  results_list <- list(
    stretch_init = stretch_init,
    stretch_lower = stretch_lower,
    stretch_upper = stretch_upper,
    shift_init = 0,
    shift_lower = shift_lower * boundary_coverage,
    shift_upper = shift_upper * boundary_coverage
  )

  return(results_list)
}

#' Optimise registration parameters with Simulated Annealing
#'
#' @param input_df TODO: Input data frame containing all replicates of gene expression in each genotype at each time point.
#' @param initial_rescale Scaling gene expression prior to registration if \code{TRUE}.
#' @param do_rescale Scaling gene expression using only overlapping time points points during registration.
#' @param min_num_overlapping_points Number of minimum overlapping time points. Shifts will be only considered if it leaves at least these many overlapping points after applying the registration function.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_ref Accession name of reference data.
#' @param num_iterations Maximum number of iterations of the algorithm. Default is 100.
#'
#' @return List of results. TODO.
#' @export
optimise_registration_params <- function(input_df, initial_guess = NA, initial_rescale = TRUE, do_rescale = FALSE,
                                         min_num_overlapping_points = 4, expression_value_threshold = 5,
                                         accession_data_to_transform, accession_data_ref,
                                         num_iterations = 100, boundary_coverage = 1) {
  # Function to optimise
  BIC_diff <- function(x, input_df, min_num_overlapping_points,
                       initial_rescale, do_rescale,
                       accession_data_to_transform, accession_data_ref) {
    stretch <- x[1]
    shift <- x[2]

    tryCatch(
      {
        BIC <- get_BIC_from_registering_data(
          input_df = input_df,
          stretches = stretch,
          shifts = shift,
          min_num_overlapping_points = min_num_overlapping_points,
          initial_rescale = initial_rescale,
          do_rescale = do_rescale,
          accession_data_to_transform = accession_data_to_transform,
          accession_data_ref = accession_data_ref,
          start_timepoint = "reference",
          optimise_shift_extreme = FALSE
        ) %>%
          suppressMessages() %>%
          suppressWarnings()

        message(stretch, " ", shift, " ", BIC)
        return(BIC)
      },
      error = function(error_message) {
        BIC <- 999

        return(BIC)
      }
    )
  }

  # Calculate boundary box and initial guess
  boundary_box <- get_boundary_box(
    input_df,
    accession_data_to_transform,
    accession_data_ref,
    min_num_overlapping_points,
    expression_value_threshold,
    boundary_coverage
  )

  if (is.na(initial_guess)) {
    stretch_init <- boundary_box$stretch_init
    shift_init <- boundary_box$shift_init
  } else {
    stretch_init <- initial_guess[1]
    shift_init <- initial_guess[2]
  }
  stretch_lower <- boundary_box$stretch_lower
  stretch_upper <- boundary_box$stretch_upper
  shift_lower <- boundary_box$shift_lower
  shift_upper <- boundary_box$shift_upper

  message("Stretch: ", round(stretch_init, 2), " in [", stretch_lower, ", ", stretch_upper, "]")
  message("Shift: ", shift_init, " in [", shift_lower, ", ", shift_upper, "]")

  # Perform SA using {optimization}
  optim_sa_res <- optimization::optim_sa(
    fun = function(x) BIC_diff(x, input_df, min_num_overlapping_points, initial_rescale, do_rescale, accession_data_to_transform, accession_data_ref),
    start = c(stretch_init, shift_init),
    trace = TRUE,
    lower = c(stretch_lower, shift_lower),
    upper = c(stretch_upper, shift_upper),
    control = list(
      t0 = 1000,
      nlimit = num_iterations,
      r = 0.85,
      rf = 3,
      dyn_rf = FALSE
    )
  )

  # Parse results
  result_df <- data.frame(
    stretch = optim_sa_res$start[1],
    shift = optim_sa_res$start[2],
    stretch_opt = round(optim_sa_res$par[1], 3),
    shift_opt = round(optim_sa_res$par[2], 3),
    value = optim_sa_res$function_value,
    is_registered = optim_sa_res$function_value < 0
  )

  results_list <- list(
    result_df = result_df,
    optim_res = optim_sa_res,
    boundary_box = boundary_box
  )

  return(results_list)
}
