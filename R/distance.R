#' Calculate distance between sample data before and after registration
#'
#' @param results Result of registration process using \code{\link{register}}.
#'
#' @return This function returns a list of data frames which includes:
#'
#' \item{registered}{distance between scaled reference and query expressions using registered time points.}
#' \item{original}{distance between scaled reference and query expressions using original time points.}
#'
#' @export
calculate_distance <- function(results) {
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

  # Retrieve data from results
  data <- results$data
  reference <- attr(data, "ref")
  query <- attr(data, "query")
  data <- unique(data[, .(expression_value = mean(expression_value)), by = .(gene_id, accession, timepoint, timepoint_reg)])

  data_query <- data[data$accession == query]
  data_ref <- data[data$accession == reference]

  # Cross join all reference and query time points
  timepoint_cj_registered <- get_timepoint_comb_registered_data(data_ref, data_query)
  timepoint_cj_original <- get_timepoint_comb_original_data(data_ref, data_query)

  # Calculate mean square distances
  dist_registered <- timepoint_cj_registered[, .(distance = mean((exp_ref - exp_query)^2)), by = .(timepoint_ref, timepoint_query)][timepoint_query >= 0]
  dist_original <- timepoint_cj_original[, .(distance = mean((exp_ref - exp_query)^2)), by = .(timepoint_ref, timepoint_query)][timepoint_query >= 0]

  # Add accession values as data attributes
  data.table::setattr(dist_registered, "ref", reference)
  data.table::setattr(dist_registered, "query", query)
  data.table::setattr(dist_original, "ref", reference)
  data.table::setattr(dist_original, "query", query)

  # Results object
  results_list <- list(
    registered = dist_registered,
    original = dist_original
  )

  return(results_list)
}

#' Cross join all original reference and query time points and expression values
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
  genes <- unique(data_query$gene_id)

  comb <- lapply(
    genes,
    function(gene) {
      comb_gene <- cross_join(
        unique(data_ref[data_ref$gene_id == gene][, .(gene_ref = gene_id, timepoint_ref = timepoint, exp_ref = expression_value)]),
        unique(data_query[data_query$gene_id == gene][, .(gene_query = gene_id, timepoint_query = timepoint, exp_query = expression_value)])
      )

      return(comb_gene)
    }
  )

  comb <- data.table::rbindlist(comb)

  # Select relevant columns
  comb <- comb[, .(gene_id = gene_ref, timepoint_ref, timepoint_query, exp_ref, exp_query)]

  return(comb)
}

#' Cross join all reference and registered query time points and expression values
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
  genes <- unique(data_query$gene_id)

  comb <- lapply(
    genes,
    function(gene) {
      comb_gene <- cross_join(
        unique(data_ref[data_ref$gene_id == gene][, .(gene_ref = gene_id, timepoint_ref = timepoint, exp_ref = expression_value)]),
        imputed_query_timepoints[imputed_query_timepoints$gene_query == gene]
      )

      return(comb_gene)
    }
  )

  comb <- data.table::rbindlist(comb)

  # Fit using cubic splines with K+3 params for each gene
  fits <- lapply(
    genes,
    function(gene) {
      fit_spline_model(data_query[data_query$gene_id == gene], x = "timepoint_reg")
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
  comb <- merge(comb, data.table::rbindlist(preds_query), by = c("gene_query", "timepoint_query"))

  # Select relevant columns
  comb <- comb[, .(gene_id = gene_ref, timepoint_ref, timepoint_query, exp_ref, exp_query)]

  return(comb)
}
