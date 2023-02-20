#' Preprocess data before registration
#'
#' \code{preprocess_data()} is a function that:
#' \item{Calculates \code{time_delta}.}
#' \item{Gets \code{mean_data}.}
#' \item{Scales data via \code{\link{scale_data}}.}
#'
#' @noRd
preprocess_data <- function(input, reference, query) {
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  accession <- NULL
  timepoint <- NULL
  expression_value <- NULL
  time_delta <- NULL

  # Make sure the data are data.tables
  all_data <- data.table::as.data.table(input)

  # Factor query and reference
  all_data[, accession := factor(accession, levels = c(reference, query), labels = c("ref", "query"))]

  # Calculate time delta for each accession
  all_data[, time_delta := timepoint - min(timepoint), by = .(accession)]

  # Get mean data
  mean_data <- data.table::copy(all_data)
  mean_data <- unique(mean_data[, .(expression_value = mean(expression_value)), by = .(gene_id, accession, timepoint, time_delta)])

  # Scale data
  scaled_data <- scale_data(mean_data, all_data)

  # Results object
  results_list <- list(
    all_data = scaled_data$all_data
  )

  return(results_list)
}

#' Scale data
#'
#' \code{scale_all_rep_data()} is a function to scale both the mean expression data and original data including all replicates.
#'
#' @param mean_data Input data containing mean of each time point.
#' @param all_data Input data including all replicates.
#'
#' @noRd
scale_data <- function(mean_data, all_data) {
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  accession <- NULL
  expression_value <- NULL
  scaled_expression_value <- NULL

  # Scale mean data
  mean_data[, scaled_expression_value := scale(expression_value, scale = TRUE, center = TRUE), by = .(gene_id, accession)]

  # Summary statistics to use for the rescaling replicates data
  gene_expression_stats <- unique(
    mean_data[, .(
      mean_val = mean(expression_value),
      sd_val = stats::sd(expression_value)
    ), by = .(gene_id, accession)]
  )

  # Left join gene_expression_stats to all_data
  all_data <- merge(
    all_data,
    gene_expression_stats,
    by = c("gene_id", "accession")
  )

  # Scale replicates data
  all_data$expression_value <- (all_data$expression_value - all_data$mean_val) / all_data$sd_val
  all_data[, c("mean_val", "sd_val") := NULL]

  # Rename expression_value for mean data
  mean_data[, c("expression_value") := NULL]
  data.table::setnames(mean_data, "scaled_expression_value", "expression_value")

  # Results object
  results_list <- list(
    all_data = all_data,
    mean_data = mean_data
  )

  return(results_list)
}
