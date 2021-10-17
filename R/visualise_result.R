#' Plot gene of interest after registration
#'
#' @param df Dataframe input after registration.
#'
#' @return Plot of gene of interest after registration.
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
      y = mean_cpm,
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
