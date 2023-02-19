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
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  accession <- NULL
  expression_value <- NULL
  registered <- NULL
  stretch <- NULL
  shift <- NULL

  # Validate parameters
  type <- match.arg(type)

  # Parse model_comparison object
  gene_facets <- results$model_comparison[, .(
    gene_id,
    gene_facet = paste0(
      gene_id, " - ", ifelse(registered, "REG", "NO_REG"),
      ifelse(
        registered,
        paste0("\n", "stretch: ", round(stretch, 2), ", shift: ", round(shift, 2)),
        ""
      )
    )
  )]

  # Left join gene_facets to data
  data <- results$data[gene_facets, on = "gene_id"]

  # Construct plot
  if (type == "registered") {
    timepoint_var <- "timepoint_reg"
    x_lab <- "Registered time"
  } else {
    timepoint_var <- "timepoint"
    x_lab <- "Time point"
  }

  gg_registered <- ggplot2::ggplot(data) +
    ggplot2::aes(
      x = !!ggplot2::sym(timepoint_var),
      y = expression_value,
      color = accession,
      fill = accession
    ) +
    ggplot2::geom_point() +
    ggplot2::stat_summary(fun = mean, geom = "line") +
    ggplot2::facet_wrap(~ gene_facet, scales = "free", ncol = ncol) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = title,
      x = x_lab,
      y = "Scaled expression"
    )

  return(gg_registered)
}
