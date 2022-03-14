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
optimise_registration_params <- function(input_df, initial_rescale = TRUE, do_rescale = FALSE,
                                         min_num_overlapping_points = 4,
                                         accession_data_to_transform, accession_data_ref,
                                         num_iterations = 100) {
  # Function to optimise
  BIC_diff <- function(x) {
    stretch <- x[1]
    shift <- x[2]

    tryCatch(
      {
        BIC <- suppressMessages(get_BIC_from_registering_data(
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
        ))

        return(BIC)
      },
      error = function(error_message) {
        BIC <- 999

        return(BIC)
      }
    )
  }

  # Stretch and shift parameters
  stretch_init <- get_approximate_stretch(
    input_df = input_df,
    accession_data_to_transform = accession_data_to_transform,
    accession_data_ref = accession_data_ref
  )
  stretch_lower <- 1
  stretch_upper <- 2 * ceiling(stretch_init) - stretch_lower

  # TODO: calculate maximum shift from data + stretch
  shift_init <- 0
  shift_lower <- -10
  shift_upper <- 10

  message("Stretch: ", stretch_init, " in [", stretch_lower, ", ", stretch_upper, "]")
  message("shift: ", shift_init, " in [", shift_lower, ", ", shift_upper, "]")

  # Perform SA using {optimization}
  optim_sa_res <- optimization::optim_sa(
    fun = BIC_diff,
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
    optim_res = optim_sa_res
  )

  return(results_list)
}
