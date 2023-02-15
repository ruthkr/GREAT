#' Summarise registration results
#'
#' @param results Registration results, output of the `register` registration process.
#'
#' @return List containing summary table, registered gene accessions, and non-registered gene accessions.
#'
#' @export
summary_registration <- function(results) {
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL

  data <- results$model_comparison

  # Summary table
  total <- nrow(data)
  reg <- sum(data$registered)
  non_reg <- total - reg

  stretch <- range(unique(data[registered == TRUE, round(stretch, 2)]))
  shift <- range(unique(data[registered == TRUE, round(shift, 2)]))

  df_summary <- data.table::data.table(
    Result = c("Total genes", "Registered genes", "Non-registered genes", "Stretch", "Shift"),
    Value = c(total,
              reg,
              non_reg,
              paste0("[", stretch[1], ", ", stretch[2], "]"),
              paste0("[", shift[1], ", ", shift[2], "]"))
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


#' Calculate prediction of reference data expression value
#'
#' @noRd
interpolate_data_ref_comparison_expression <- function(data_to_transform_time, data_ref_dt) {
  # Suppress "no visible binding for global variable" note
  shifted_time <- NULL

  # Calculate diff
  data_ref_dt$diff <- data_ref_dt$shifted_time - data_to_transform_time

  # If outside of comparable range (time is smaller than all data_ref_dt time or bigger than all)
  if (all(data_ref_dt$diff > 0) | all(data_ref_dt$diff < 0)) {
    return(NA)
  }

  # Otherwise, cut down data_ref observations to the two nearest timepoints to the data_to_transform time
  data_ref_dt$diff <- abs(data_ref_dt$shifted_time - data_to_transform_time)
  data.table::setorder(data_ref_dt, diff)
  nearest_points <- data_ref_dt[1:2, ]

  # Linearly interpolate between these points to estimate the comparison expression value
  data.table::setorder(nearest_points, shifted_time) # so [1] is earlier time
  time_diff <- nearest_points$shifted_time[2] - nearest_points$shifted_time[1]
  expression_diff <- nearest_points$expression_value[2] - nearest_points$expression_value[1]
  grad <- expression_diff / time_diff
  pred_expression <- nearest_points$expression_value[1] + (nearest_points$diff[1]) * grad

  return(pred_expression)
}
