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

#' Visualise distances between samples from different time points
#'
#' @description
#' Function `plot_heatmap()` allows users to plot distances between samples from
#' different time points to investigate the similarity of progression of gene
#' expression states between species before or after registration.
#'
#' @param data Input data frame containing sample distance between two different species.
#' @param match_timepoints If \code{TRUE}, will match query time points to reference time points.
#' @param title Optional plot title.
#' @param axis_fontsize Font size of X and Y axes labels.
#'
#' @return Distance heatmap of gene expression profiles over time between two different species.
#'
#' @export
plot_heatmap <- function(data, match_timepoints = FALSE, title = NULL, axis_fontsize = NULL) {
  # Suppress "no visible binding for global variable" note
  timepoint_ref <- NULL
  timepoint_query <- NULL
  distance <- NULL

  # Retrieve accession values from results
  reference <- attr(data, "ref")
  query <- attr(data, "query")

  # Synchronise time points
  if (match_timepoints) {
    equal_timepoints <- intersect(data$timepoint_ref, data$timepoint_query)
    data <- data[timepoint_query %in% equal_timepoints & timepoint_ref %in% equal_timepoints]
  }

  # Construct plot
  gg_distance <- ggplot2::ggplot(data) +
    ggplot2::aes(
      x = as.factor(timepoint_query),
      y = as.factor(timepoint_ref),
      fill = log(distance)
    ) +
    ggplot2::geom_tile() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = axis_fontsize),
      axis.text.y = ggplot2::element_text(size = axis_fontsize),
      panel.border = ggplot2::element_blank(),
      legend.position = "top",
      legend.justification = "right",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, -10, -10),
      legend.title = ggplot2::element_text(size = 10),
      legend.key.height = ggplot2::unit(5, "pt")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(label.position = "top")
    ) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(
      title = title,
      x = query,
      y = reference
    )

  # Synchronise time points
  if (match_timepoints) {
    gg_distance <- gg_distance + ggplot2::coord_fixed()
  }

  return(gg_distance)
}
