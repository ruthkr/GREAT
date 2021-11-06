#' Summarise registration results
#'
#' @param model_comparison Input data frame, element \code{model_comparison} of result list of \code{scale_and_register_data()}.
#'
#' @return List containing summary table, registered gene accessions, and non-registered gene accessions.
#' @importFrom rlang .data
#' @export
summary_model_comparison <- function(model_comparison) {
  # Suppress "no visible binding for global variable" note
  gene <- NULL

  # Summary table
  total <- nrow(model_comparison)
  reg <- sum(model_comparison$BIC_registered_is_better)
  non_reg <- total - reg

  stretch <- model_comparison %>%
    dplyr::filter(.data$BIC_registered_is_better == TRUE) %>%
    dplyr::pull(.data$stretch) %>%
    unique() %>%
    sort()

  shift <- model_comparison %>%
    dplyr::filter(.data$BIC_registered_is_better == TRUE) %>%
    dplyr::pull(.data$shift) %>%
    unique() %>%
    sort()

  df_summary <- data.frame(
    Result = c("Total genes", "Registered genes", "Non-registered genes", "Stretch", "Shift"),
    Value = c(total, reg, non_reg, paste(stretch, collapse = ", "), paste0("[", min(shift), ", ", max(shift), "]"))
  )

  # List of registered and non-registered genes
  registered_genes <- unique(model_comparison[model_comparison$BIC_registered_is_better, gene])
  non_registered_genes <- unique(model_comparison[!model_comparison$BIC_registered_is_better, gene])

  # Results object
  results_list <- list(
    df_summary = df_summary,
    registered_genes = registered_genes,
    non_registered_genes = non_registered_genes
  )

  return(results_list)
}
