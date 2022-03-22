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

  stretch_lower <- 0.5 * stretch_init
  # stretch_upper <- 1.5 * ceiling(stretch_init) - stretch_lower
  stretch_upper <- 1.5 * ceiling(stretch_init)

  # Shift limits
  all_data_df <- data.table::as.data.table(input_df)

  mean_df <- get_mean_data(
    exp = all_data_df,
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
  ) %>%
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

#' Optimise registration parameters with Simulated Annealing for single gene
#' @noRd
#' @importFrom rlang .data
optimise_registration_params_single_gene <- function(input_df,
                                                     initial_guess = NA,
                                                     initial_rescale = FALSE,
                                                     do_rescale = TRUE,
                                                     min_num_overlapping_points = 4,
                                                     expression_value_threshold = 5,
                                                     accession_data_to_transform,
                                                     accession_data_ref,
                                                     start_timepoint,
                                                     is_data_normalised,
                                                     num_iterations = 100,
                                                     boundary_coverage = 1) {
  # Function to optimise
  BIC_diff <- function(x) {
    stretch <- x[1]
    shift <- x[2]

    tryCatch(
      {
        BIC <- get_BIC_from_registering_data(
          input_df = input_df,
          stretches = stretch,
          shifts = shift,
          min_num_overlapping_points = min_num_overlapping_points,
          expression_value_threshold = expression_value_threshold,
          initial_rescale = initial_rescale,
          do_rescale = do_rescale,
          accession_data_to_transform = accession_data_to_transform,
          accession_data_ref = accession_data_ref,
          start_timepoint = start_timepoint,
          is_data_normalised = is_data_normalised,
          optimise_shift_extreme = FALSE
        ) %>%
          suppressMessages() %>%
          suppressWarnings()

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
  ) %>%
    suppressMessages()

  if (any(is.na(initial_guess))) {
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

  cli::cli_alert_info("Stretch: { round(stretch_init, 2) } in [{ stretch_lower }, { stretch_upper }]")
  cli::cli_alert_info("Shift: { round(shift_init, 2) } in [{ shift_lower }, { shift_upper }]")

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
  locus_name <- unique(input_df$locus_name)

  result_df <- data.frame(
    gene = locus_name,
    stretch = round(optim_sa_res$par[1], 3),
    shift = round(optim_sa_res$par[2], 3),
    BIC_diff = optim_sa_res$function_value,
    is_registered = optim_sa_res$function_value < 0
  )

  trace_df <- optim_sa_res$trace %>%
    as.data.frame() %>%
    dplyr::mutate(
      gene = locus_name,
      is_registered = .data$loss < 0
    ) %>%
    dplyr::filter(.data$is_registered) %>%
    dplyr::select(
      .data$gene,
      stretch = .data$x_1,
      shift = .data$x_2,
      BIC_diff = .data$loss,
      .data$is_registered
    ) %>%
    dplyr::distinct()

  results_list <- list(
    optimum_params_df = result_df,
    candidate_params_df = trace_df
  )

  return(results_list)
}

#' Optimise registration parameters with Simulated Annealing
#'
#' @param input_df Input data frame containing all replicates of gene expression in each genotype at each time point.
#' @param genes List of genes to optimise.
#' @param initial_rescale Scaling gene expression prior to registration if \code{TRUE}.
#' @param do_rescale Scaling gene expression using only overlapping time points points during registration.
#' @param min_num_overlapping_points Number of minimum overlapping time points. Shifts will be only considered if it leaves at least these many overlapping points after applying the registration function.
#' @param expression_value_threshold Expression value threshold. Remove expressions if maximum is less than the threshold. If \code{NULL} keep all data.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_ref Accession name of reference data.
#' @param num_iterations Maximum number of iterations of the algorithm. Default is 100.
#' @param boundary_coverage Coverage factor of the boundary box. Default is 1.
#'
#' @return List of optimum registration parameters, \code{optimum_params_df}, and other candidate registration parameters, \code{candidate_params_df} for all genes.
#' @export
#' @importFrom rlang .data
optimise_registration_params <- function(input_df,
                                         genes = NULL,
                                         initial_rescale = FALSE,
                                         do_rescale = TRUE,
                                         min_num_overlapping_points = 4,
                                         expression_value_threshold = 5,
                                         accession_data_to_transform,
                                         accession_data_ref,
                                         start_timepoint,
                                         is_data_normalised,
                                         num_iterations = 100,
                                         boundary_coverage = 1) {
  # Validate genes
  if (is.null(genes)) {
    genes <- unique(input_df$locus_name)
  }

  # Apply optimise_registration_params_single_gene() over all genes
  raw_results <- genes %>%
    purrr::map(
      function(gene) {
        curr_df <- input_df %>%
          dplyr::filter(.data$locus_name == gene)

        opt_res <- optimise_registration_params_single_gene(
          input_df = curr_df,
          initial_guess = NA,
          initial_rescale,
          do_rescale,
          min_num_overlapping_points,
          expression_value_threshold,
          accession_data_to_transform,
          accession_data_ref,
          start_timepoint,
          is_data_normalised,
          num_iterations,
          boundary_coverage
        )

        return(opt_res)
      }
    )

  # Parse raw results
  optimum_params_reduced <- raw_results %>%
    purrr::map(~ purrr::pluck(.x, "optimum_params_df")) %>%
    purrr::reduce(dplyr::bind_rows)

  candidate_params_reduced <- raw_results %>%
    purrr::map(~ purrr::pluck(.x, "candidate_params_df")) %>%
    purrr::reduce(dplyr::bind_rows)

  results_list <- list(
    optimum_params_df = optimum_params_reduced,
    candidate_params_df = candidate_params_reduced
  )

  return(results_list)
}
