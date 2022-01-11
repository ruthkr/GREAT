#' Plot gene of interest after registration
#'
#' @param reg_result_df Data frame of registration results, output from registration process.
#' @param model_comparison_df Data frame of model comparison, also output from registration process.
#' @param gene_accession List of gene accessions, default is \code{first_genes} which will take first 25 genes.
#' @param title Optional plot title.
#' @param ncol Number of columns in the plot grid. By default this is calculated automatically.
#' @param sync_timepoints Whether to synchronise maximum time points for each accession, by default \code{FALSE}.
#'
#' @return Plot of gene of interest after registration process.
#' @importFrom rlang .data
#' @export
plot_registration_results <- function(reg_result_df, model_comparison_df = NULL, gene_accession = "first_genes", title = NULL, ncol = NULL, sync_timepoints = FALSE) {
  # Make sure that the accession is in character format
  reg_result_df$accession <- as.character(reg_result_df$accession)

  # Filter gene using given gene of interests
  if (gene_accession == "first_genes") {
    first_25_genes <- reg_result_df %>%
      dplyr::pull(.data$locus_name) %>%
      unique() %>%
      utils::head(25)

    reg_result_df <- reg_result_df %>%
      dplyr::filter(.data$locus_name %in% first_25_genes)
  } else {
    reg_result_df <- reg_result_df %>%
      dplyr::filter(.data$locus_name %in% gene_accession)
  }

  # Synchronise maximum time points for each accession
  # TODO: consider n additional time points
  if (sync_timepoints) {
    max_timepoints <- reg_result_df %>%
      dplyr::filter(
        !is.na(.data$shifted_time),
        !is.na(.data$expression_value)
      ) %>%
      dplyr::group_by(.data$locus_name, .data$accession) %>%
      dplyr::summarise(
        max_timepoint = max(.data$shifted_time),
        .groups = "drop"
      ) %>%
      dplyr::group_by(.data$locus_name) %>%
      dplyr::summarise(
        max_timepoint = min(.data$max_timepoint),
        .groups = "drop"
      )

    reg_result_df <- dplyr::left_join(
      reg_result_df,
      max_timepoints,
      by = "locus_name"
    ) %>%
      dplyr::filter(.data$shifted_time <= .data$max_timepoint) %>%
      dplyr::select(-.data$max_timepoint)
  }

  if (!is.null(model_comparison_df)) {
    reg_result_df <- reg_result_df %>%
      dplyr::left_join(
        model_comparison_df %>%
          dplyr::select(locus_name = .data$gene, .data$stretch, .data$shift),
        by = "locus_name"
      ) %>%
      dplyr::mutate(
        is_registered = ifelse(.data$is_registered, "REG", "NO_REG"),
        locus_name = paste0(
          .data$locus_name, " - ", .data$is_registered, "\n",
          "stretch: ", round(.data$stretch, 2),
          ", shift: ", round(.data$shift, 2)
        )
      )
  }

  # Plot
  gg_registered <- ggplot2::ggplot(reg_result_df) +
    ggplot2::aes(
      x = .data$shifted_time,
      y = .data$expression_value,
      color = .data$accession,
      fill = .data$accession,
      # TODO: group = interaction(locus_name, bra_gene)
    ) +
    ggplot2::geom_point(size = 0.4) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ .data$locus_name, scales = "free", ncol = ncol) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = title,
      x = "Registered time (d)",
      y = "Normalised expression"
    )

  # TODO: handle replicates
  # if (FALSE) {
  #   gg_registered <- gg_registered +
  #     ggplot2::stat_summary(fun = mean, geom = "line", size = 1) +
  #     ggplot2::stat_summary(
  #       fun.data = mean_se,
  #       fun.args = list(mult = 1),
  #       geom = "ribbon",
  #       color = NA,
  #       alpha = 0.3
  #     )
  # }

  return(gg_registered)
}

#' Visualise distances between samples from different time points
#'
#' @description
#' Function `plot_heatmap()` allows users to plot distances between samples from different time points to investigate the similarity of progression of gene expression states between species before or after registration.
#'
#' @param sample_dist_df Input data frame contains sample distance between two different species.
#' @param title Optional plot title.
#' @param axis_fontsize Font size of X and Y axes labels.
#' @param same_min_timepoint If \code{FALSE}, the default, will not take data with the same minimum time point.
#' @param same_max_timepoint If \code{FALSE}, the default, will not take data with the same maximum time point.
#'
#' @return Distance heatmap of gene expression profiles over time between two different species.
#' @importFrom rlang .data
#' @export
plot_heatmap <- function(sample_dist_df, title = NULL, axis_fontsize = NULL, same_min_timepoint = FALSE, same_max_timepoint = FALSE) {

  # Define timepoint x and y
  sample_dist_df <- sample_dist_df %>%
    dplyr::mutate(
      timepoint_x = as.numeric(stringr::str_extract(.data$x_sample, "(?<=-)\\d+")),
      timepoint_y = as.numeric(stringr::str_extract(.data$y_sample, "(?<=-)\\d+"))
    )

  # Take data with the same minimum time points
  if (same_min_timepoint) {
    sample_dist_df <- sample_dist_df %>%
      dplyr::filter(
        .data$timepoint_x >= min(.data$timepoint_y)
      )
  }

  # Take data with the same maximum time points
  if (same_max_timepoint) {
    sample_dist_df <- sample_dist_df %>%
      dplyr::filter(
        .data$timepoint_x <= max(.data$timepoint_y)
      )
  }

  # Change class of x_sample and y_sample as factor
  sample_dist_df$x_sample <- factor(sample_dist_df$x_sample, levels = unique(sort(sample_dist_df$x_sample)))
  sample_dist_df$y_sample <- factor(sample_dist_df$y_sample, levels = unique(sort(sample_dist_df$y_sample)))


  p <- ggplot2::ggplot(sample_dist_df) +
    ggplot2::aes(
      x = .data$x_sample,
      y = .data$y_sample,
      fill = log(.data$distance)
    ) +
    ggplot2::geom_tile() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, size = axis_fontsize),
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
    viridis::scale_fill_viridis() +
    ggplot2::labs(
      title = title,
      x = "",
      y = ""
    )

  return(p)
}
