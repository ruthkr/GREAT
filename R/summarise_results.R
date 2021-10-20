#' Summarise registration results
#'
#' @param model_comparison Input data frame, element \code{model_comparison} of result list of \code{scale_and_register_data()}.
#'
#' @return List containing summary table, registered gene accessions, and non-registered gene accessions.
#' @export
summary_model_comparison <- function(model_comparison) {
  # Summary table
  total <- nrow(model_comparison)
  reg <- length(model_comparison$BIC_registered_is_better == TRUE)
  non_reg <- total - reg

  df_summary <- data.frame(
    Result = c("Total genes", "Registered genes", "Non-registered genes"),
    Count = c(total, reg, non_reg)
  )

  # List of registered and non-registered genes
  registered_genes <- unique(model_comparison[model_comparison$BIC_registered_is_better, gene])
  non_registered_genes <- unique(model_comparison[!model_comparison$BIC_registered_is_better, gene])

  # Results object
  results_list <- list(
    summary = df_summary,
    registered_genes = registered_genes,
    non_registered_genes = non_registered_genes
  )

  return(results_list)
}
