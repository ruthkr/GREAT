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

  # Calculate BIC for Hypothesis 2
  BIC_separate <- compare_dynamics_H2(all_data)

  # Apply registration
  if (!optimise_registration_parameters) {
    # Check that stretches and shifts are numeric
    if (any(is.na(stretches)) | any(is.na(shifts))) {
      stop(cli::format_error(c(
        "{.var stretches} and {.var shifts} must be numeric vectors",
        "x" = "You supplied vectors with {.cls NA} values."
      )))
    }

    # Apply registration
    mean_data_reg <- apply_registration(mean_data, stretches, shifts)
    all_data_reg <- apply_registration(all_data, stretches, shifts)

    BIC_combined <- compare_dynamics_H1(all_data_reg)
    message("Registered: ", BIC_separate > BIC_combined)
  }


  # Final data processing
  # TODO: undo factor() of the reference and query for interpretability
  # TODO: (?) accession counts summary with {cli}

  results_list <- list(
    all_data = data.table::data.table(all_data),
    mean_data = data.table::data.table(mean_data),
    all_data_reg = data.table::data.table(all_data_reg),
    mean_data_reg = data.table::data.table(mean_data_reg)
  )

  return(results_list)
}
