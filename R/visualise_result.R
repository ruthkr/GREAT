#' Plot gene of interest after registration
#'
#' @param df Dataframe input after registration.
#'
#' @return Plot of gene of interest after registration.
#' @export
plot_registered_gene_of_interest <- function(df,
                                             gene_accession = "all") {

  # Make sure that the accession is in character format
  df$accession <- as.character(df$accession)

  # Filter gene using given gene of interests
  if (gene_accession == "all"){
    df <- df
  } else {
    df <- df %>%
      dplyr::filter(accession %in% gene_accession)
  }

  # Plot
  gg_registered <- ggplot2::ggplot(df) +
    ggplot2::aes(
      x = shifted_time,
      y = expression_value,
      color = accession,
      fill = accession
    ) +
    # ggplot2::stat_summary(fun = mean, geom = "line", size = 1) +
    # ggplot2::stat_summary(
    #   fun.data = mean_se,
    #   fun.args = list(mult = 1),
    #   geom = "ribbon",
    #   color = NA,
    #   alpha = 0.3
    # ) +
    ggplot2::geom_point(size = 0.4) +
    ggplot2::geom_line() +
    ggplot2::xlab("Registered time (d)") +
    ggplot2::ylab("Normalised expression") +
    ggplot2::facet_wrap(~locus_name, scales = "free", ncol = 2) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 10),
      axis.text = ggplot2::element_text(size = 6),
      strip.text = ggplot2::element_text(face = "italic"),
      legend.margin = ggplot2::margin(21, 0, 0, 0)
    )

  return(gg_registered)
}

#' Visualise distances between samples from different time points to investigate the similarity of progression of gene expression states between species before or after registration
#'
#' @param sample_dist_df Input data frame contains sample distance between two different species.
#' @param ylabel Label on y axis.
#' @param y_axis_fontsize Font size of y axis label.
#' @param same_min_timepoint
#'
#' @return Distance heatmap of gene expression profiles over time between two different species.
#' @export
make_heatmap <- function(sample_dist_df, ylabel, y_axis_fontsize = 6, same_min_timepoint = TRUE) {

  if (same_min_timepoint) {
    sample_dist_df <- sample_dist_df %>%
      dplyr::mutate(
        timepoint_x = stringr::str_extract(x_sample, "(?<=-)\\d+") %>%
          as.numeric(),
        timepoint_y = stringr::str_extract(y_sample, "(?<=-)\\d+") %>%
          as.numeric()
      ) %>%
      dplyr::filter(timepoint_x >= min(timepoint_y))
  } else {
    sample_dist_df <- sample_dist_df
  }

  # Change class of x_sample and y_sample as factor
  sample_dist_df$x_sample <- factor(sample_dist_df$x_sample, levels = unique(sort(sample_dist_df$x_sample)))
  sample_dist_df$y_sample <- factor(sample_dist_df$y_sample, levels = unique(sort(sample_dist_df$y_sample)))

  p <- ggplot2::ggplot(sample_dist_df) +
    ggplot2::aes(
      x = x_sample,
      y = y_sample,
      fill = log(distance)
    ) +
    ggplot2::geom_tile() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = 6
      ),
      axis.text.y = ggplot2::element_text(size = y_axis_fontsize),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
      plot.margin = ggplot2::margin(0, 0, -10, 0),
      panel.background = ggplot2::element_blank(),
      legend.position = "top",
      legend.justification = "right",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, -10, -10),
      legend.text = ggplot2::element_text(size = 4, vjust = -0.5),
      legend.title = ggplot2::element_text(size = 8),
      legend.key.height = ggplot2::unit(0.2, "cm"),
    ) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(label.position = "top")) +
    viridis::scale_fill_viridis() +
    ggplot2::ylab(ylabel) +
    ggplot2::xlab("")

  return(p)
}


