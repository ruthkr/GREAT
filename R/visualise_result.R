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

#' Make data heatmaps
make_data_heatmaps <- function(D.mean, D.scaled, D.registered, D.scaled.NR, D.scaled.R, D.registered.R) {
  p.mean <- make_heatmap(D.mean, ylabel = "mean expression")
  p.scaled <- make_heatmap(D.scaled, ylabel = "scaled mean expression")
  p.registered <- make_heatmap(D.registered, ylabel = "registered & scaled mean expression")

  p.scaled.NR <- make_heatmap(D.scaled.onlyNR, ylabel = "scaled mean expression (only not registered genes)")
  p.scaled.R <- make_heatmap(D.scaled.onlyR, ylabel = "scaled mean expression (only registered genes)")
  p.registered.R <- make_heatmap(D.registered.R, ylabel = "registered & scaled mean expression (only registered genes)")

  p.all <- cowplot::plot_grid(p.mean, p.scaled.NR, p.scaled, p.scaled.R, p.registered, p.registered.R, ncol = 2)

  return(p.all)
}

#' Make heatmap
make_heatmap <- function(D, ylabel, y.axis.fontsize = 6) {
  D$x.sample <- factor(D$x.sample, levels = unique(sort(D$x.sample)))
  D$y.sample <- factor(D$y.sample, levels = unique(sort(D$y.sample)))

  p <- ggplot2::ggplot(D) +
    ggplot2::aes(
      x = x.sample,
      y = y.sample,
      fill = log(distance)
    ) +
    ggplot2::geom_tile() +
    # scale_fill_viridis()+
    # theme_classic()+
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = 6
      ),
      axis.text.y = ggplot2::element_text(size = y.axis.fontsize),
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
    ggplot2::facet_wrap(~title, nrow = 1) +
    ggplot2::ylab(ylabel) +
    ggplot2::xlab("")
  # ggplot2::ggtitle(title)

  return(p)
}


#' Make all heatmaps
make_heatmap_all <- function(D, title) {
  D$x.sample <- factor(D$x.sample, levels = unique(sort(D$x.sample)))
  D$y.sample <- factor(D$y.sample, levels = unique(sort(D$y.sample)))

  p <- ggplot2::ggplot(D) +
    ggplot2::aes(
      x = x.sample,
      y = y.sample,
      fill = log(distance)
    ) +
    ggplot2::geom_tile() +
    # viridis::scale_fill_viridis()+
    ggplot2::theme_classic() +
    ggplot2::facet_wrap(~title, ncol = 1, scales = "free") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = 6
      ),
      axis.text.y = ggplot2::element_text(size = 6),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
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
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::ggtitle(title)

  return(p)
}

