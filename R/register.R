register <- function(input,
                     stretches = NA,
                     shifts = NA,
                     reference,
                     query,
                     overlapping_percent = 0.5,
                     optimise_registration_parameters = FALSE,
                     optimisation_config = list(num_iterations = 60)) {
  # Validate column names
  match_names(
    x = colnames(input),
    lookup = c("gene_id", "accession", "timepoint", "expression_value", "replicate"),
    error = "Must review the column names of your input data:",
    name_string = "column names"
  )

  # Validate accession names
  match_names(
    x = unique(input$accession),
    lookup = c(reference, query),
    error = "Must review the supplied {.var reference} and {.var query} values:",
    name_string = "accession values"
  )

  # Preprocess data
  processed_data <- preprocess_data(input, reference, query)

  all_data <- processed_data$all_data
  mean_data <- processed_data$mean_data

  # TODO: match_starting_timepoints()
  # @return data.table shifted data, delta_timepoint_ref, delta_timepoint_query

  # TODO: apply_registration() for a single gene

  # Final data processing
  # TODO: undo factor() of the reference and query for interpretability
  # TODO: (?) accession counts summary with {cli}

  results_list <- list(
    all_data = all_data,
    mean_data = mean_data
  )

  return(results_list)
}
