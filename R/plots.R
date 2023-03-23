#' Plot gene of interest after registration
#'
#' @param results Registration results, output of the \code{\link{register}} registration process.
#' @param type Type of plot, determines whether to use "registered" or "original" time points. By default, "registered".
#' @param genes_list Optional vector indicating the \code{gene_id} values to be plotted.
#' @param title Optional plot title.
#' @param ncol Number of columns in the plot grid. By default this is calculated automatically.
#'
#' @return Plot of genes of interest after registration process (\code{type = "registered"}) or showing original time points (\code{type = "original"}).
#'
#' @export
plot_registration_results <- function(results,
                                      type = c("registered", "original"),
                                      genes_list = NULL,
                                      title = NULL,
                                      ncol = NULL) {
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  accession <- NULL
  timepoint_reg <- NULL
  expression_value <- NULL

  # Validate parameters
  type <- match.arg(type)

  # Parse results
  data <- results$data
  model_comparison <- results$model_comparison
  reference <- attr(data, "ref")
  query <- attr(data, "query")

  # Select genes to be plotted
  genes <- unique(data[, gene_id])

  if (any(!is.null(genes_list))) {
    if (!inherits(genes_list, "character")) {
      stop(
        cli::format_error(c(
          "{.var genes_list} must be a {.cls character} vector.",
          "x" = "You supplied vectors with {.cls {class(genes_list)}} values."
        )),
        call. = FALSE
      )
    }

    data <- data[data$gene_id %in% genes_list]
    model_comparison <- model_comparison[model_comparison$gene_id %in% genes_list]
  } else if (length(genes) > 50) {
    cli::cli_alert_info("The first 25 genes will be shown. To override this, use the {.var genes_list} parameter.")
    model_comparison <- model_comparison[model_comparison$gene_id %in% genes[1:25]]
  }

  # Parse model_comparison object
  gene_facets <- parse_gene_facets(model_comparison, type)
  data <- data[gene_facets, on = "gene_id"]

  # Plot labels
  if (type == "registered") {
    timepoint_var <- "timepoint_reg"
    x_lab <- "Registered time"
  } else {
    timepoint_var <- "timepoint"
    x_lab <- "Time point"
  }

  # Construct plot
  gg_registered <- ggplot2::ggplot(data) +
    ggplot2::aes(
      x = !!ggplot2::sym(timepoint_var),
      y = expression_value,
      color = accession
      # fill = accession
    ) +
    ggplot2::geom_point() +
    # ggplot2::stat_summary(fun = mean, geom = "line") +
    ggplot2::facet_wrap(~gene_facet, scales = "free", ncol = ncol) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = title,
      x = x_lab,
      y = "Scaled expression"
    )

  # Add model curve layers
  if (type == "registered") {
    # Count registered and unregistered genes
    registered_count <- length(model_comparison[model_comparison$registered, gene_id])
    unregistered_count <- length(model_comparison[!model_comparison$registered, gene_id])

    # Get curves for H1
    if (registered_count > 0) {
      preds_H1 <- get_H1_model_curves(data, model_comparison)
      preds_H1 <- merge(preds_H1, gene_facets, by = "gene_id")
    }
    # Get curves for H2
    if (unregistered_count > 0) {
      preds_H2 <- get_H2_model_curves(data, model_comparison, reference, query)
      preds_H2 <- merge(preds_H2, gene_facets, by = "gene_id")
    }

    # Bind predictions
    if (unregistered_count == 0) {
      preds <- preds_H1
    } else if (registered_count == 0) {
      preds <- preds_H2
    } else {
      preds <- rbind(preds_H1, preds_H2)
    }

    gg_registered <- gg_registered +
      ggplot2::geom_line(
        mapping = ggplot2::aes(
          x = timepoint_reg,
          y = expression_value,
          group = interaction(accession, gene_id)
        ),
        data = preds,
        linetype = "dashed",
        linewidth = 0.5
      )
  }

  # Fix legend order
  gg_registered <- gg_registered +
    ggplot2::scale_color_manual(
      breaks = c(reference, query),
      values = scales::hue_pal()(2)
    )

  return(gg_registered)
}

#' Parse \code{gene_facet} faceting variable for plotting
#'
#' @noRd
parse_gene_facets <- function(model_comparison, type) {
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  registered <- NULL
  stretch <- NULL
  shift <- NULL

  if (type == "registered") {
    gene_facets <- model_comparison[, .(
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
  } else {
    gene_facets <- model_comparison[, .(
      gene_id,
      gene_facet = gene_id
    )]
  }

  return(gene_facets)
}

#' Get curves for Hypothesis H1
#'
#' @noRd
get_H1_model_curves <- function(data, model_comparison) {
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  accession <- NULL
  timepoint_reg <- NULL

  # Get registered genes only
  genes <- unique(model_comparison[model_comparison$registered, gene_id])

  # Fit using cubic splines with K+3 params for each gene
  fits <- lapply(
    genes,
    function(gene) {
      fit_spline_model(
        data[data$gene_id == gene],
        x = "timepoint_reg"
      )
    }
  )
  names(fits) <- genes

  # Predict query expression values
  preds <- data.table::rbindlist(
    lapply(
      genes,
      function(gene) {
        data <- unique(data[data$gene_id == gene][, .(gene_id, timepoint_reg)])
        data <- data.table::data.table(
          gene_id = gene,
          timepoint_reg = seq(min(data$timepoint_reg), max(data$timepoint_reg), 1)
        )
        data[, .(gene_id, timepoint_reg, expression_value = stats::predict(fits[gene][[1]], newdata = data))]
      }
    )
  )

  preds[, accession := character()][]

  return(preds)
}

#' Get curves for Hypothesis H2
#'
#' @noRd
get_H2_model_curves <- function(data, model_comparison, reference, query) {
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  accession <- NULL
  timepoint_reg <- NULL

  # Get unregistered genes only
  genes <- unique(model_comparison[!model_comparison$registered, gene_id])

  # Get reference and query data
  data_ref <- data[data$accession == reference]
  data_query <- data[data$accession == query]

  # Fit using cubic splines with K+3 params for each gene
  fits_ref <- lapply(
    genes,
    function(gene) {
      fit_spline_model(
        data_ref[data_ref$gene_id == gene],
        x = "timepoint_reg"
      )
    }
  )
  fits_query <- lapply(
    genes,
    function(gene) {
      fit_spline_model(
        data_query[data_query$gene_id == gene],
        x = "timepoint_reg"
      )
    }
  )
  names(fits_ref) <- genes
  names(fits_query) <- genes

  # Predict query expression values
  preds_ref <- data.table::rbindlist(
    lapply(
      genes,
      function(gene) {
        data <- unique(data_ref[data_ref$gene_id == gene][, .(gene_id, timepoint_reg)])
        data <- data.table::data.table(
          gene_id = gene,
          timepoint_reg = seq(min(data$timepoint_reg), max(data$timepoint_reg), 1)
        )
        data[, .(gene_id, timepoint_reg, expression_value = stats::predict(fits_ref[gene][[1]], newdata = data))]
      }
    )
  )

  preds_query <- data.table::rbindlist(
    lapply(
      genes,
      function(gene) {
        data <- unique(data_query[data_query$gene_id == gene][, .(gene_id, timepoint_reg)])
        data <- data.table::data.table(
          gene_id = gene,
          timepoint_reg = seq(min(data$timepoint_reg), max(data$timepoint_reg), 1)
        )
        data[, .(gene_id, timepoint_reg, expression_value = stats::predict(fits_query[gene][[1]], newdata = data))]
      }
    )
  )

  # Combine reference and query curves
  preds <- rbind(
    preds_ref[, accession := reference],
    preds_query[, accession := query]
  )[]

  return(preds)
}

#' Visualise distances between samples from different time points
#'
#' \code{plot_heatmap()} is a function that allows users to plot distances
#' between samples from different time points to investigate the similarity of
#' progression of gene expression states between species before or after
#' registration.
#'
#' @param results Results containing distances between two different reference and query data, output of \code{\link{calculate_distance}}.
#' @param type Type of plot, determines whether to use "registered" or "original" time points. By default, "registered".
#' @param match_timepoints If \code{TRUE}, will match query time points to reference time points.
#' @param title Optional plot title.
#' @param axis_fontsize Font size of X and Y axes labels.
#'
#' @return Distance heatmap of gene expression profiles over time between two different species.
#'
#' @export
plot_heatmap <- function(results,
                         type = c("registered", "original"),
                         match_timepoints = FALSE,
                         title = NULL,
                         axis_fontsize = NULL) {
  # Suppress "no visible binding for global variable" note
  timepoint_ref <- NULL
  timepoint_query <- NULL
  distance <- NULL

  # Validate parameters
  type <- match.arg(type)
  if (type == "registered") {
    data <- results$registered
  } else {
    data <- results$original
  }

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
