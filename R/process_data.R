#' Preprocess data before registration
#'
#' \code{preprocess_data()} is a function that:
#' \item{Calculates \code{time_delta}.}
#' \item{Calculates expression \code{var} values for each timepoint.}
#' \item{Scales data via \code{\link{scale_data}}.}
#'
#' @noRd
preprocess_data <- function(input, reference, query, exp_sd = NA, scaling_method = c("none", "z-score", "min-max")) {
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

  # Filter genes that are missing one accession
  all_data <- filter_incomplete_accession_pairs(all_data)

  # Filter genes that do not change expression over time (sd == 0)
  all_data <- filter_unchanged_expressions(all_data)

  cli::cli_alert_info("Will process {length(unique(all_data$gene_id))} gene{?s}.")

  # Calculate time delta for each accession
  all_data[, time_delta := timepoint - min(timepoint), by = .(gene_id, accession)]

  # Calculate expression variance
  all_data <- calc_variance(all_data, exp_sd)

  # Scale data
  scaled_data <- scale_data(all_data, scaling_method)

  return(scaled_data)
}

#' Filter genes that are missing one accession
#'
#' \code{filter_incomplete_accession_pairs()} is a function to filter out genes that are missing one accession.
#'
#' @param all_data Input data including all replicates.
#'
#' @noRd
filter_incomplete_accession_pairs <- function(all_data) {
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  count <- NULL

  # Detect genes with only one accession
  accession_counts <- unique(all_data, by = c("gene_id", "accession"))[, .(count = .N), by = .(gene_id)]
  discarded_genes <- accession_counts[accession_counts$count == 1]$gene_id

  n_genes <- length(discarded_genes)
  if (n_genes > 0) {
    cli::cli_alert_warning("{n_genes} gene{?s} only {?has/have} data with one accession, therefore {?it/they} will be discarded.")
  }

  return(all_data[!gene_id %in% discarded_genes])
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
#' @param all_data Input data including all replicates.
#' @param scaling_method Scaling method applied to data prior to registration process. Either \code{none} (default), \code{z-score}, or \code{min-max}.
#'
#' @noRd
scale_data <- function(all_data, scaling_method = c("none", "z-score", "min-max")) {
  # Validate parameters
  scaling_method <- match.arg(scaling_method)

  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  accession <- NULL
  expression_value <- NULL
  scaled_expression_value <- NULL

  if (scaling_method == "z-score") {
    # Calculate mean and standard deviation of expression in all_data by accession
    all_data[,
      c("mean_val", "sd_val") := .(mean(expression_value), stats::sd(expression_value)),
      by = .(gene_id, accession)
    ]

    # Scale replicates data
    all_data$expression_value <- (all_data$expression_value - all_data$mean_val) / all_data$sd_val
    all_data$var <- pmax(all_data$var / (all_data$sd_val)^2, 0.75^2)
    all_data[, c("mean_val", "sd_val") := NULL]
  } else if (scaling_method == "min-max") {
    # Calculate minimum and maximum of expression in all_data by accession
    all_data[,
      c("min_val", "max_val") := .(min(expression_value), max(expression_value)),
      by = .(gene_id, accession)
    ]

    # Scale replicates data
    all_data$expression_value <- (all_data$expression_value - all_data$min_val) / (all_data$max_val - all_data$min_val)
    all_data$var <- all_data$var / (all_data$max_val - all_data$min_val)^2
    all_data[, c("min_val", "max_val") := NULL]
  }

  return(all_data)
}
