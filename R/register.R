#' Register or synchronize different expression profiles
#'
#' `register()` is a function to register expression profiles a user
#' wishes to compare.
#'
#' @param input Input data frame containing all replicates of gene expression in each genotype at each time point.
#' @param stretches Candidate registration stretch factors to apply to query data.
#' @param shifts Candidate registration shift values to apply to query data.
#' @param reference Accession name of reference data.
#' @param query Accession name of query data.
#' @param overlapping_percent Number of minimum overlapping time points. Shifts will be only considered if it leaves at least these many overlapping points after applying the registration function.
#' @param optimise_registration_parameters Whether to optimise registration parameters with Simulated Annealing. By default, \code{FALSE}.
#' @param optimisation_config List with optional arguments to modify the Simulated Annealing optimisation.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load a data frame from the sample data
#' all_data <- system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "greatR") %>%
#'   utils::read.csv()
#'
#' # Running the registration
#' registration_results <- register(
#'   input = all_data,
#'   reference = "Ro18",
#'   query = "Col0",
#'   stretches = 2.3,
#'   shifts = 1.6
#' )
#' }
register <- function(input,
                     stretches = NA,
                     shifts = NA,
                     reference,
                     query,
                     overlapping_percent = 0.5,
                     optimise_registration_parameters = TRUE,
                     optimisation_config = list(num_iterations = 60)) {
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  accession <- NULL
  timepoint <- NULL
  timepoint_reg <- NULL
  expression_value <- NULL
  replicate <- NULL

  # Validate column names
  match_names(
    x = colnames(input),
    lookup = c("gene_id", "accession", "timepoint", "expression_value", "replicate"),
    error = "Must review the column names of your input data:",
    name_string = "column names"
  )

  # Validate accession names
  match_names(
    x = c(reference, query),
    lookup = unique(input$accession),
    error = "Must review the supplied {.var reference} and {.var query} values:",
    name_string = "accession values"
  )

  # Preprocess data
  processed_data <- preprocess_data(input, reference, query)
  all_data <- processed_data$all_data
  gene_id_list <- unique(all_data$gene_id)
  cli::cli_alert_info("Will process {length(gene_id_list)} gene{?s}.")

  # Begin registration logic
  if (optimise_registration_parameters) {
    # Registration with optimisation
    cli::cli_h1("Starting registation with optimisation")
    cli::cli_alert_info("Using computed stretches and shifts search space limits. User-defined parameters will be ignored.")

    # Run optimisation
    results <- lapply(
      cli::cli_progress_along(
        x = gene_id_list,
        format = "{cli::pb_spin} Optimising registration parameters for genes ({cli::pb_current}/{cli::pb_total}) [{cli::pb_elapsed_clock}]",
        format_done = "{cli::col_green(cli::symbol$tick)} Optimising registration parameters for genes ({cli::pb_total}/{cli::pb_total}) {cli::col_white(paste0('[', cli::pb_elapsed, ']'))}",
        clear = FALSE
      ),
      function(gene_index) {
        # Filter single gene data
        gene_data <- all_data[all_data$gene_id == gene_id_list[gene_index]]

        # Calculate BIC for Hypothesis 2
        loglik_separate <- calc_loglik_H2(gene_data)

        # Register for Hypothesis 1
        results <- register_with_optimisation(gene_data, loglik_separate, overlapping_percent, optimisation_config)

        return(results)
      }
    )
  } else {
    cli::cli_h1("Starting manual registration")
    # Check that stretches and shifts are numeric
    if (any(is.na(stretches), is.na(shifts))) {
      stop(cli::format_error(c(
        "{.var stretches} and {.var shifts} must be numeric vectors",
        "x" = "You supplied vectors with {.cls NA} values."
      )))
    }

    # Apply manual registration
    results <- lapply(
      cli::cli_progress_along(
        x = gene_id_list,
        format = "{cli::pb_spin} Applying registration for genes ({cli::pb_current}/{cli::pb_total}) [{cli::pb_elapsed_clock}]",
        format_done = "{cli::col_green(cli::symbol$tick)} Applying registration for genes ({cli::pb_total}/{cli::pb_total}) {cli::col_white(paste0('[', cli::pb_elapsed, ']'))}",
        clear = FALSE
      ),
      function(gene_index) {
        # Filter single gene data
        gene_data <- all_data[all_data$gene_id == gene_id_list[gene_index]]

        # Calculate BIC for Hypothesis 2
        loglik_separate <- calc_loglik_H2(gene_data)

        # Register for Hypothesis 1
        results <- register_manually(gene_data, stretches, shifts, loglik_separate)
        # TODO: generalise for multiple stretches and shifts
        # results <- get_best_parameters()

        return(results)
      }
    )
  }

  # Aggregate results
  all_data_reg <- Reduce(rbind, lapply(results, function(x) x$data_reg))
  model_comparison <- Reduce(rbind, lapply(results, function(x) x$model_comparison))

  # Left join registered time points
  setnames(all_data_reg, "timepoint", "timepoint_reg")
  all_data <- merge(
    all_data,
    all_data_reg,
    by = c("gene_id", "accession", "expression_value", "replicate")
  )

  # Restore original query and reference accession names
  all_data[, c("time_delta") := NULL]
  all_data[, accession := factor(accession, levels = c("ref", "query"), labels = c(reference, query))][, .(gene_id, accession, timepoint, timepoint_reg, expression_value, replicate)]

  # Add accession values as data attributes
  data.table::setattr(all_data, "ref", reference)
  data.table::setattr(all_data, "query", query)

  # Results object
  results_list <- list(
    data = all_data,
    model_comparison = model_comparison
  )

  return(results_list)
}

#' @noRd
register_with_optimisation <- function(data, loglik_separate, overlapping_percent, optimisation_config) {
  # Run optimisation
  optimised_params <- optimise(data, overlapping_percent, optimisation_config)

  # Apply registration
  stretches <- optimised_params$par[1]
  shifts <- optimised_params$par[2]
  data_reg <- apply_registration(data, stretches, shifts)

  # Calculate model comparison
  loglik_combined <- optimised_params$function_value

  # Model comparison
  model_comparison <- compare_H1_and_H2(data_reg, stretches, shifts, loglik_combined, loglik_separate)

  # Results object
  results_list <- list(
    data_reg = data_reg,
    model_comparison = model_comparison
  )

  return(results_list)
}

#' @noRd
register_manually <- function(data, stretch, shift, loglik_separate) {
  # Apply registration
  data_reg <- apply_registration(data, stretch, shift)

  # Calculate model comparison
  loglik_combined <- calc_loglik_H1(data_reg)

  # Model comparison
  model_comparison <- compare_H1_and_H2(data_reg, stretch, shift, loglik_combined, loglik_separate)

  # Results object
  results_list <- list(
    data_reg = data_reg,
    model_comparison = model_comparison
  )

  return(results_list)
}
