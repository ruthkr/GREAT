#' Calculate distance between sample data before and after registration
#'
#' @param registration_results Result of registration process using \code{\link{scale_and_register_data}}.
#' @param gene_col Column name of gene accession, default is \code{locus_name}.
#' @param compare_ref_vs_transform If \code{TRUE}, the default, only comparison between reference data and data to transform is considered.
#' @param accession_data_ref Accession name of reference data.
#'
#' @return This function returns a list of data frames which includes:
#' \item{distance_mean_df}{distance of mean expression values.}
#' \item{distance_scaled_mean_df}{distance of scaled mean expression (all genes).}
#' \item{distance_scaled_mean_df_only_nonreg}{distance of scaled mean expression (only non-registered genes).}
#' \item{distance_scaled_mean_df_only_reg}{distance of scaled mean expression (only registered genes).}
#' \item{distance_registered_df}{distance of registered & scaled mean expression (all genes).}
#' \item{distance_registered_df_only_reg}{distance of registered & scaled mean expression (only registered genes).}
#'
#' @export
calculate_between_sample_distance <- function(registration_results,
                                              gene_col = "locus_name",
                                              compare_ref_vs_transform = TRUE,
                                              accession_data_ref) {
  # Parse registration_results list
  mean_df <- registration_results$mean_df
  mean_df_sc <- registration_results$mean_df_sc
  imputed_mean_df <- registration_results$imputed_mean_df

  # Convert all to wide format ready for distance calculation
  # mean_df
  mean.dt.w <- reformat_for_distance_calculation(
    dt = mean_df,
    sample_id_cols = c("accession", "timepoint"),
    gene_col = gene_col,
    expression_col = "expression_value"
  )

  # normalised mean_df
  mean.dt.sc.w <- reformat_for_distance_calculation(
    dt = mean_df_sc,
    sample_id_cols = c("accession", "timepoint"),
    gene_col = gene_col,
    expression_col = "sc.expression_value"
  )

  # imputed_mean_df - all genes
  imputed.mean.dt.w <- reformat_for_distance_calculation(
    imputed_mean_df,
    sample_id_cols = c("accession", "shifted_time"),
    gene_col = gene_col,
    expression_col = "expression_value"
  )

  # TODO: change `condition == FALSE/TRUE` by `isFALSE(condition)/isTRUE(condition)` or `condition/!condition`
  # same, but for subsets of REGISTERED / NOT REGISTERED genes.
  # distance between samples, only using genes which are found best model is not registered
  not.registered.genes <- unique(imputed_mean_df$locus_name[imputed_mean_df$is_registered == FALSE])
  mean.dt.sc.w.not.registered <- mean.dt.sc.w[mean.dt.sc.w$locus_name %in% not.registered.genes, ]

  # distance between samples, only using genes which are found best when ARE registered
  registered.genes <- unique(imputed_mean_df$locus_name[imputed_mean_df$is_registered == TRUE])
  mean.dt.sc.w.registered <- mean.dt.sc.w[mean.dt.sc.w$locus_name %in% registered.genes, ]

  # after registration, but only for registered genes
  imputed.mean.dt.w.registered <- imputed.mean.dt.w[imputed.mean.dt.w$locus_name %in% registered.genes, ]

  # Calculate distance between each sample in each of the wide format tables.
  # shifting genes results in NAs for some genes at some timepoints.
  # therefore different numbers of dimensions (genes) depending on the comparison.
  # therefore euclidean distance is not appropriate.
  # Distance used is mean of [squared distance between each gene / absolute mean of expression of that gene in the sample]. Is calculated for each gene for which have data in both samples.
  # Mean of these values (divided by number of genes is calculated for) is reported.
  D.mean <- calc_sample_distance(
    mean.dt.w,
    gene_col = gene_col,
    compare_ref_vs_transform = compare_ref_vs_transform,
    accession_data_ref = accession_data_ref
  )

  D.scaled <- calc_sample_distance(
    mean.dt.sc.w,
    gene_col = gene_col,
    compare_ref_vs_transform = compare_ref_vs_transform,
    accession_data_ref = accession_data_ref
  )

  D.registered <- calc_sample_distance(
    imputed.mean.dt.w,
    gene_col = gene_col,
    compare_ref_vs_transform = compare_ref_vs_transform,
    accession_data_ref = accession_data_ref
  )

  D.scaled.not.registered.genes <- calc_sample_distance(
    mean.dt.sc.w.not.registered,
    gene_col = gene_col,
    compare_ref_vs_transform = compare_ref_vs_transform,
    accession_data_ref = accession_data_ref
  )

  D.scaled.registered.genes <- calc_sample_distance(
    mean.dt.sc.w.registered,
    gene_col = gene_col,
    compare_ref_vs_transform = compare_ref_vs_transform,
    accession_data_ref = accession_data_ref
  )

  D.registered.registered.genes <- calc_sample_distance(
    imputed.mean.dt.w.registered,
    gene_col = gene_col,
    compare_ref_vs_transform = compare_ref_vs_transform,
    accession_data_ref = accession_data_ref
  )

  # Titles for heatmaps with shared scales
  attr(D.mean, "title") <- "mean expression"
  attr(D.scaled, "title") <- "scaled mean expression (all genes)"
  attr(D.registered, "title") <- "registered & scaled mean expression (all genes)"
  attr(D.scaled.not.registered.genes, "title") <- "scaled mean expression (only not-registered genes)"
  attr(D.scaled.registered.genes, "title") <- "scaled mean expression (only registered genes)"
  attr(D.registered.registered.genes, "title") <- "registered & scaled mean expression (only registered genes)"

  # Results object
  results_list <- list(
    "distance_mean_df" = D.mean,
    "distance_scaled_mean_df" = D.scaled,
    "distance_scaled_mean_df_only_nonreg" = D.scaled.not.registered.genes,
    "distance_scaled_mean_df_only_reg" = D.scaled.registered.genes,
    "distance_registered_df" = D.registered,
    "distance_registered_df_only_reg" = D.registered.registered.genes
  )

  return(results_list)
}

#' Reformat distance calculation
#'
#' @noRd
reformat_for_distance_calculation <- function(dt, sample_id_cols, gene_col, expression_col) {
  # Avoid dcast() "Aggregate function missing, defaulting to 'length'" error
  dt <- unique(dt)

  # Concatenate sample.id columns to generate sample ids
  dt$sample.id <- dt[[sample_id_cols[1]]]
  if (length(sample_id_cols) > 1) {
    for (i in 2:length(sample_id_cols)) {
      dt$sample.id <- paste0(dt[["sample.id"]], "-", dt[[sample_id_cols[i]]])
    }
  }

  # Subset to just the relevant columns
  dt <- subset(dt, select = c("sample.id", gene_col, expression_col))

  # Convert to wide format
  dt.w <- data.table::dcast(dt, locus_name ~ sample.id, value.var = eval(expression_col))

  # Remove columns that contain only NAs
  dt.w <- dt.w[, colSums(is.na(dt.w)) != nrow(dt.w), with = FALSE]

  return(dt.w)
}

#' Calculate sample distance wrapper
#'
#' @noRd
calc_sample_distance <- function(df, gene_col, compare_ref_vs_transform = TRUE, accession_data_ref) {
  data.cols <- names(df)[names(df) != eval(gene_col)]

  # Allocate size for results
  size <- length(data.cols)
  i.cols <- numeric(length = size * (size - 1))
  j.cols <- numeric(length = size * (size - 1))
  ds <- numeric(length = size * (size - 1))
  count <- 1

  for (i in 1:(size - 1)) {
    i.col <- data.cols[i]
    for (j in (i + 1):size) {
      j.col <- data.cols[j]

      curr.dt <- subset(df, select = c(i.col, j.col))
      dist <- calculate_pairwise_sample_distance_main(curr.dt)

      # Distance metric is symmetrical
      i.cols[2 * count - 1] <- j.cols[2 * count] <- i.col
      j.cols[2 * count - 1] <- i.cols[2 * count] <- j.col
      ds[c(2 * count - 1, 2 * count)] <- dist
      count <- count + 1
    }
  }

  if (compare_ref_vs_transform) {
    out_df <- data.table::data.table(
      data.frame(
        x_sample = i.cols,
        y_sample = j.cols,
        distance = ds
      )
    ) %>%
      dplyr::filter(
        stringr::str_extract(.data$y_sample, "^.*?(?=-)") != stringr::str_extract(.data$x_sample, "^.*?(?=-)"),
        stringr::str_extract(.data$y_sample, "^.*?(?=-)") == accession_data_ref
      )
  } else {
    out_df <- data.table::data.table(
      data.frame(
        x_sample = i.cols,
        y_sample = j.cols,
        distance = ds
      )
    )
  }

  return(out_df)
}

#' Calculate pairwise sample distance
#'
#' @noRd
calculate_pairwise_sample_distance_main <- function(df) {
  # Omit NA in the data
  df <- stats::na.omit(df)

  # Calculate distance
  df$sq_diff <- (df[, 1] - df[, 2])^2

  d <- mean(df$sq_diff)

  return(d)
}
