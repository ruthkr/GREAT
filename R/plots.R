#' Plot gene of interest after registration
#'
#' @param results Registration results, output of the `register` registration process.
#' @param type Type of plot, determines whether to use "registered" or "original" time points. By default, "registered".
#' @param title Optional plot title.
#' @param ncol Number of columns in the plot grid. By default this is calculated automatically.
#'
#' @return Plot of genes of interest after registration process.
#'
#' @export
plot_registration_results <- function(results, type = c("registered", "original"), title = NULL, ncol = NULL) {
  # Validate parameters
  type <- match.arg(type)

  # Parse model_comparison object
  # reg_result_df <- reg_result_df %>%
  #   dplyr::left_join(
  #     model_comparison_df %>%
  #       dplyr::select(locus_name = .data$gene, .data$stretch, .data$shift),
  #     by = "locus_name"
  #   ) %>%
  #   dplyr::mutate(
  #     is_registered = ifelse(.data$is_registered, "REG", "NO_REG"),
  #     locus_name = paste0(
  #       .data$locus_name, " - ", .data$is_registered, "\n",
  #       "stretch: ", round(.data$stretch, 2),
  #       ", shift: ", round(.data$shift, 2)
  #     )
  #   )

  # Construct plot
  if (type == "registered") {
    timepoint_var <- "timepoint_reg"
    x_lab <- "Registered time"
  } else {
    timepoint_var <- "timepoint"
    x_lab <- "Time point"
  }

  gg_registered <- ggplot2::ggplot(results$data) +
    ggplot2::aes_string(
      x = timepoint_var,
      y = "expression_value",
      color = "accession",
      fill = "accession"
    ) +
    ggplot2::geom_point() +
    ggplot2::stat_summary(fun = mean, geom = "line") +
    ggplot2::facet_wrap(~gene_id, scales = "free", ncol = ncol) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = title,
      x = x_lab,
      y = "Scaled expression"
    )

  return(gg_registered)
}
