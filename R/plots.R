#' Visualise registration results
#'
#' @param x Input object.
#'  - For [plot.res_greatR()]: registration results, output of the [register()] registration process.
#'  - For [plot.dist_greatR()]: pairwise distances between reference and query time points, output of [calculate_distance()].
#' @param type Type of plot, whether to use registration "result" or "original" time points. By default, "result".
#' @param genes_list Optional vector indicating the \code{gene_id} values to be plotted.
#' @param show_rep_mean Whether to show \code{replicate} mean values.
#' @param match_timepoints If \code{TRUE}, will match query time points to reference time points.
#' @param ncol Number of columns in the plot grid. By default this is calculated automatically.
#' @param title Optional plot title.
#' @param ... Arguments to be passed to methods (ignored).
#'
#' @name plot
#' @return
#'  - For [plot.res_greatR()]: plot of genes of interest after registration process (\code{type = "result"}) or showing original time points (\code{type = "original"}).
#'  - For [plot.dist_greatR()]: distance heatmap of gene expression profiles over time between reference and query.
NULL
# > NULL

#' @rdname plot
#' @export
plot.res_greatR <- function(x,
                            type = c("result", "original"),
                            genes_list = NULL,
                            show_rep_mean = FALSE,
                            ncol = NULL,
                            title = NULL,
                            ...) {
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  accession <- NULL
  timepoint_reg <- NULL
  expression_value <- NULL

  # Validate parameters
  type <- match.arg(type)

  # Parse results
  data <- x$data
  model_comparison <- x$model_comparison
  reference <- attr(data, "ref")
  query <- attr(data, "query")
  scaling_method <- x$fun_args$scaling_method

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
  if (type == "result") {
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
    ) +
    ggplot2::geom_point() +
    {
      if (show_rep_mean) ggplot2::stat_summary(fun = mean, geom = "line")
    } +
    ggplot2::facet_wrap(~gene_facet, scales = "free", ncol = ncol) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
    theme_greatR() +
    ggplot2::labs(
      title = title,
      x = x_lab,
      y = ifelse(scaling_method == "none", "Expression", "Scaled expression"),
      colour = NULL
    )

  # Add model curve layers
  if (type == "result") {
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
      values = greatR_palettes$disc
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
  BIC_diff <- NULL

  if (type == "result") {
    gene_facets <- model_comparison[, .(
      gene_id,
      gene_facet = paste0(
        gene_id, " - ", ifelse(registered, "REG", "NON-REG"),
        ifelse(
          registered,
          paste0("\n", "BIC diff: ", round(BIC_diff, 2), ", stretch: ", round(stretch, 2), ", shift: ", round(shift, 2)),
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
          timepoint_reg = seq(min(data$timepoint_reg), max(data$timepoint_reg), length.out = 25)
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
          timepoint_reg = seq(min(data$timepoint_reg), max(data$timepoint_reg), length.out = 25)
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
          timepoint_reg = seq(min(data$timepoint_reg), max(data$timepoint_reg), length.out = 25)
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

#' @rdname plot
#' @export
plot.dist_greatR <- function(x,
                             type = c("result", "original"),
                             match_timepoints = TRUE,
                             title = NULL,
                             ...) {
  # Suppress "no visible binding for global variable" note
  timepoint_ref <- NULL
  timepoint_query <- NULL
  distance <- NULL

  # Validate parameters
  type <- match.arg(type)
  if (type == "result") {
    data <- x$result
  } else {
    data <- x$original
  }

  # Retrieve accession values from results
  reference <- attr(data, "ref")
  query <- attr(data, "query")

  # Synchronise time points
  if (all(match_timepoints, type == "result")) {
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
    theme_greatR() +
    ggplot2::theme(
      legend.position = "top",
      legend.justification = "right",
      legend.margin = ggplot2::margin(0, 0, -5, 0),
      legend.title = ggplot2::element_text(size = 10),
      legend.key.height = ggplot2::unit(5, "pt")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(label.position = "top")
    ) +
    ggplot2::scale_fill_gradientn(colours = greatR_palettes$cont) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(
      title = title,
      x = query,
      y = reference
    )

  # Synchronise time points
  if (all(match_timepoints, type == "result")) {
    gg_distance <- gg_distance + ggplot2::coord_fixed()
  }

  return(gg_distance)
}

#' @noRd
theme_greatR <- function(base_size = 12, base_family = "sans") {
  (
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(
        size = ceiling(base_size * 0.7),
        colour = "black"
      ),
      axis.title = ggplot2::element_text(
        size = ceiling(base_size * 0.8)
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        fill = NA,
        colour = "black",
        linewidth = 1,
        linetype = "solid"
      ),
      strip.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(
        size = ceiling(base_size * 0.75)
      ),
      legend.position = "right",
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      plot.title = ggplot2::element_text(
        size = ceiling(base_size * 1.1), face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(
        size = ceiling(base_size * 1.05)
      )
    )
  )
}

#' @noRd
greatR_palettes <- list(
  cont = c("#3e1190", "#4268b7", "#4195bd", "#51bdb9", "#8fdc9f", "#f9ef93"),
  disc = c("#1b9e77", "#f38400"),
  hist = c("REG" = "#4268b7", "NON-REG" = "#d98484")
)
