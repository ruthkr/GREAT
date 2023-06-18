#' Preprocess data before registration
#'
#' \code{preprocess_data()} is a function that:
#' \item{Calculates \code{time_delta}.}
#' \item{Gets \code{mean_data}.}
#' \item{Scales data via \code{\link{scale_data}}.}
#'
#' @noRd
preprocess_data <- function(input, reference, query, scaling_method = c("scale", "normalise")) {
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  accession <- NULL
  timepoint <- NULL
  expression_value <- NULL
  time_delta <- NULL

  # Validate parameters
  scaling_method <- match.arg(scaling_method)

  # Make sure the data are data.tables
  all_data <- data.table::as.data.table(input)

  # Factor query and reference
  all_data[, accession := factor(accession, levels = c(reference, query), labels = c("ref", "query"))]

  # Filter genes that do not change expression over time (sd == 0)
  all_data <- filter_unchanged_expressions(all_data)

  # Calculate time delta for each accession
  all_data[, time_delta := timepoint - min(timepoint), by = .(gene_id, accession)]

  # Get mean data
  mean_data <- data.table::copy(all_data)
  mean_data <- unique(mean_data[, .(expression_value = mean(expression_value)), by = .(gene_id, accession, timepoint, time_delta)])

  # Scale data
  scaled_data <- scale_data(mean_data, all_data, scaling_method)

  # Results object
  results_list <- list(
    all_data = scaled_data$all_data
  )

  return(results_list)
}

#' Filter genes that do not change over time
#'
#' \code{filter_unchanged_expressions()} is a function to filter out genes that do not change expression over time.
#'
#' @param all_data Input data including all replicates.
#'
#' @noRd
filter_unchanged_expressions <- function(all_data) {
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  accession <- NULL
  expression_value <- NULL
  sd <- NULL
  exp_sd <- NULL

  # Calculate standard deviations
  gene_sds <- all_data[, .(exp_sd = sd(expression_value)), by = .(gene_id, accession)]
  discarded_genes <- unique(gene_sds[exp_sd <= 1e-16]$gene_id)
  n_genes <- length(discarded_genes)
  if (n_genes > 0) {
    cli::cli_alert_warning("{n_genes} gene{?s} {?does/do} not change over time, therefore {?it/they} will be discarded.")
  }

  return(all_data[!gene_id %in% discarded_genes])
}

#' Scale data
#'
#' \code{scale_all_rep_data()} is a function to scale both the mean expression data and original data including all replicates.
#'
#' @param mean_data Input data containing mean of each time point.
#' @param all_data Input data including all replicates.
#' @param scaling_method Scaling method applied to data prior to registration process. Either \code{scale} (default), or \code{normalise}.
#'
#' @noRd
scale_data <- function(mean_data, all_data, scaling_method = c("scale", "normalise")) {
  # Validate parameters
  scaling_method <- match.arg(scaling_method)

  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  accession <- NULL
  expression_value <- NULL
  scaled_expression_value <- NULL

  if (scaling_method == "scale") {
    # Scale mean data
    mean_data[, scaled_expression_value := scale(expression_value, scale = TRUE, center = FALSE), by = .(gene_id, accession)]

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
    all_data$expression_value <- all_data$expression_value / all_data$sd_val
    all_data[, c("mean_val", "sd_val") := NULL]
  } else if (scaling_method == "normalise") {
    # Summary statistics to use for the rescaling replicates data
    gene_expression_stats <- mean_data[, .(
       min_expression_value = min(expression_value),
       max_expression_value = max(expression_value)
     ),
     by = .(gene_id, accession)
    ]

    # Left join gene_expression_stats to all_data
    all_data <- merge(
      all_data,
      gene_expression_stats,
      by = c("gene_id", "accession")
    )

    # Scale replicates data
    all_data$expression_value <- (all_data$expression_value - all_data$min_expression_value) / (all_data$max_expression_value - all_data$min_expression_value)
    all_data[, c("min_expression_value", "max_expression_value") := NULL]
  }

  # Results object
  results_list <- list(
    all_data = all_data
  )

  return(results_list)
}
