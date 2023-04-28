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

  # Filter out not differentially expressed genes (sd == 0)
  all_data <- get_diff_expressed_data(all_data)

  # Calculate time delta for each accession
  all_data[, time_delta := timepoint - min(timepoint), by = .(gene_id, accession)]

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

#' Get differentially expressed data
#'
#' \code{get_diff_expressed_data()} is a function to filter out not differentially expressed genes.
#'
#' @param all_data Input data including all replicates.
#'
#' @noRd
get_diff_expressed_data <- function(all_data) {
  gene_sds <- all_data[, .(exp_sd = sd(expression_value)), by = .(gene_id, accession)]
  discarded_genes <- unique(gene_sds[exp_sd <= 1e-16]$gene_id)
  n_genes <- length(discarded_genes)
  if (n_genes > 0) {
    cli::cli_alert_warning("{n_genes} gene{?s} {?is/are} not differentially expressed, therefore discarded.")
  }

  return(all_data[!gene_id %in% discarded_genes])
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

  scaling_method <- "original"
  # scaling_method <- "ratio_gene"
  # scaling_method <- "ratio_gene_accession"
  # scaling_method <- "normalise"


  if (scaling_method == "original") {
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
    # all_data$expression_value <- (all_data$expression_value - all_data$mean_val) / all_data$sd_val
    all_data$expression_value <- all_data$expression_value / all_data$sd_val
    all_data[, c("mean_val", "sd_val") := NULL]
  } else if (scaling_method == "ratio_gene") {
    # Scale mean data
    mean_data[, scaled_expression_value := scale(expression_value, scale = TRUE, center = FALSE), by = .(gene_id)]

    # Summary statistics to use for the rescaling replicates data
    gene_expression_stats <- mean_data[, .(
      scaling = mean(scaled_expression_value / expression_value, na.rm = TRUE)
    ),
    by = .(gene_id)
    ]

    # Left join gene_expression_stats to all_data
    all_data <- merge(
      all_data,
      gene_expression_stats,
      by = c("gene_id")
    )

    # Scale replicates data
    all_data$expression_value <- (all_data$expression_value) * all_data$scaling
    all_data[, "scaling" := NULL]
  } else if (scaling_method == "ratio_gene_accession") {
    # Scale mean data
    mean_data[, scaled_expression_value := scale(expression_value, scale = TRUE, center = FALSE), by = .(gene_id, accession)]

    # Summary statistics to use for the rescaling replicates data
    gene_expression_stats <- mean_data[, .(
      scaling = mean(scaled_expression_value / expression_value, na.rm = TRUE)
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
    all_data$expression_value <- (all_data$expression_value) * all_data$scaling
    all_data[, "scaling" := NULL]
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

  # Rename expression_value for mean data
  # mean_data[, c("expression_value") := NULL]
  # data.table::setnames(mean_data, "scaled_expression_value", "expression_value")

  # Results object
  results_list <- list(
    all_data = all_data
    # mean_data = mean_data
  )

  return(results_list)
}
