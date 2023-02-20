#' Calculate distance between sample data before and after registration
#'
#' @param results Result of registration process using \code{\link{register}}.
#' @param type Type of comparison, determines whether to use "registered" or "original" time points. By default, "registered".
#'
#' @return This function returns a list of data frames which includes:
#' \item{distance_mean_df}{distance of mean expression values.}
#' \item{distance_scaled_mean_df}{distance of scaled mean expression (all genes).}
#' \item{distance_scaled_mean_df_only_nonreg}{distance of scaled mean expression (only non-registered genes).}
#' \item{distance_scaled_mean_df_only_reg}{distance of scaled mean expression (only registered genes).}
#' \item{distance_registered_df}{distance of registered & scaled mean expression (all genes).}
#' \item{distance_registered_df_only_reg}{distance of registered & scaled mean expression (only registered genes).}
#'
#' @export
calculate_distance <- function(results,
                                              type = c("registered", "original")) {
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  accession <- NULL
  timepoint <- NULL
  timepoint_reg <- NULL
  timepoint_ref <- NULL
  timepoint_query <- NULL
  expression_value <- NULL
  exp_ref <- NULL
  exp_query <- NULL

  # Validate parameters
  type <- match.arg(type)

  # Retrieve data from results
  data <- results$data
  reference <- attr(data, "ref")
  query <- attr(data, "query")
  data <- unique(data[, .(expression_value = mean(expression_value)), by = .(gene_id, accession, timepoint, timepoint_reg)])

  data_query <- data[data$accession == query]
  data_ref <- data[data$accession == reference]

  if (type == "registered") {
    timepoint_cross_join <- get_timepoint_comb_registered_data(data_ref, data_query)
  } else {
    timepoint_cross_join <- get_timepoint_comb_original_data(data_ref, data_query)
  }

  # Calculate mean square distances
  distances <- timepoint_cross_join[, .(distance = mean((exp_ref - exp_query)^2)), by = .(timepoint_ref, timepoint_query)]

  return(distances)
}

#' Cross join all original reference and query timepoints and expression values
#'
#' @noRd
get_timepoint_comb_original_data <- function(data_ref, data_query) {
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  gene_ref <- NULL
  gene_query <- NULL
  timepoint <- NULL
  timepoint_ref <- NULL
  timepoint_query <- NULL
  expression_value <- NULL
  exp_ref <- NULL
  exp_query <- NULL

  # Perform cross join
  comb <- cross_join(
    unique(data_ref[, .(gene_ref = gene_id, timepoint_ref = timepoint, exp_ref = expression_value)]),
    unique(data_query[, .(gene_query = gene_id, timepoint_query = timepoint, exp_query = expression_value)])
  )[gene_ref == gene_query]

  # Select relevant columns
  comb <- comb[, .(gene_id = gene_ref, timepoint_ref, timepoint_query, exp_ref, exp_query)]

  return(comb)
}

#' Cross join all reference and registered query timepoints and expression values
#'
#' @noRd
get_timepoint_comb_registered_data <- function(data_ref, data_query) {
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  gene_ref <- NULL
  gene_query <- NULL
  timepoint <- NULL
  timepoint_reg <- NULL
  timepoint_ref <- NULL
  timepoint_query <- NULL
  expression_value <- NULL
  exp_ref <- NULL
  exp_query <- NULL

  # The imputed query time points to estimate expression values for
  timepoint_ranges_query <- data_query[, .(min_t = ceiling(min(timepoint_reg)), max_t = floor(max(timepoint_reg))), by = "gene_id"]

  imputed_query_timepoints <- data.table::rbindlist(
    Map(
      function(x, min_t, max_t) {
        data.table::data.table(gene_query = rep(x, max_t - min_t + 1), timepoint_query = seq(min_t, max_t))
      }, timepoint_ranges_query$gene_id, timepoint_ranges_query$min_t, timepoint_ranges_query$max_t
    )
  )

  # Perform cross join
  comb <- cross_join(
    unique(data_ref[, .(gene_ref = gene_id, timepoint_ref = timepoint, exp_ref = expression_value)]),
    imputed_query_timepoints
  )[gene_ref == gene_query]

  # Fit using cubic splines with K+3 params for each gene
  genes <- unique(data_query$gene_id)
  fits <- lapply(
    genes,
    function(gene) {
      stats::lm(
        expression_value ~ splines::bs(timepoint_reg, df = 6, degree = 3),
        data = data_query[data_query$gene_id == gene]
      )
    }
  )
  names(fits) <- genes

  # Predict query expression values
  preds_query <- lapply(
    genes,
    function(gene) {
      data <- unique(comb[comb$gene_query == gene][, .(timepoint_reg = timepoint_query)])
      data[, .(gene_query = gene, timepoint_query = timepoint_reg, exp_query = stats::predict(fits[gene][[1]], newdata = data))]
    }
  )

  # Left join to cross join
  comb <- merge(comb, Reduce(rbind, preds_query), by = c("gene_query", "timepoint_query"))

  # Select relevant columns
  comb <- comb[, .(gene_id = gene_ref, timepoint_ref, timepoint_query, exp_ref, exp_query)]

  return(comb)
}
