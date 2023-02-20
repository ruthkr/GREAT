#' Summarise registration results
#'
#' @param results Registration results, output of the \code{\link{register}} registration process.
#'
#' @return List containing summary table, registered gene accessions, and non-registered gene accessions.
#'
#' @export
summary_registration <- function(results) {
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL

  # Summarise results
  data <- results$model_comparison

  total <- nrow(data)
  reg <- sum(data$registered)
  non_reg <- total - reg

  stretch <- range(unique(data[data$registered, round(stretch, 2)]))
  shift <- range(unique(data[data$registered, round(shift, 2)]))
  stretch_range <- paste0("[", stretch[1], ", ", stretch[2], "]")
  shift_range <- paste0("[", shift[1], ", ", shift[2], "]")

  # Create summary table
  df_summary <- data.table::data.table(
    Result = c("Total genes", "Registered genes", "Non-registered genes", "Stretch", "Shift"),
    Value = c(total, reg, non_reg, stretch_range, shift_range)
  )

  # List of registered and non-registered genes
  registered_genes <- unique(data[data$registered, gene_id])
  non_registered_genes <- unique(data[!data$registered, gene_id])

  # Results object
  results_list <- list(
    df_summary = df_summary,
    registered_genes = registered_genes,
    non_registered_genes = non_registered_genes
  )

  return(results_list)
}
