#' Calculate mean expression values from all expression data with replicates
#'
#' @param exp Input data frame contains all replicates of gene expression in each genotype at each timepoint.
#' @param expression_value_threshold Expression value threshold. Remove expressions if maximum is less than the threshold. If \code{NULL} keep all data.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param is_data_normalised TRUE if dataset has been normalised prior to registration process.
#'
#' @return A data frame contains only mean expression data.
#' @export
get_mean_data <- function(exp,
                          expression_value_threshold = 5,
                          accession_data_to_transform,
                          is_data_normalised = FALSE) {
  # Suppress "no visible binding for global variable" note
  mean_expression_value <- NULL
  expression_value <- NULL
  locus_name <- NULL
  accession <- NULL
  tissue <- NULL
  timepoint <- NULL
  keep <- NULL
  keep_final <- NULL

  # Calculate mean of each timepoint by adding a column called "expression_value"
  # TODO: make vector in mean_df a non-hardcoded parameter
  exp[, mean_expression_value := mean(expression_value), by = list(locus_name, accession, tissue, timepoint)]

  mean_df <- unique(exp[, c("locus_name", "accession", "tissue", "timepoint", "mean_expression_value")])

  data_ref_df <- mean_df[mean_df$accession != accession_data_to_transform]

  if (!is_data_normalised) {
    # Filter mean_df to remove genes with expression lower than the threshold, and less than half timepoints expressed greater than 1
    if (!is.null(expression_value_threshold)) {
      data_ref_df[, keep := (max(mean_expression_value) > expression_value_threshold | mean(mean_expression_value > 1) > 0.5), by = .(locus_name)]
    } else {
      data_ref_df[, keep := mean(mean_expression_value) > 0, by = .(locus_name)]
    }
  } else {
    data_ref_df[, keep := TRUE]
  }

  keep_data_ref_genes <- unique(data_ref_df$locus_name[data_ref_df$keep == TRUE])
  discard_data_ref_genes <- unique(data_ref_df$locus_name[data_ref_df$keep == FALSE])

  # Filter mean_df to remove all data to transform genes with all zeros values
  data_to_transform_df <- mean_df[mean_df$locus_name %in% keep_data_ref_genes & mean_df$accession == accession_data_to_transform]
  data_to_transform_df[, keep_final := (mean(mean_expression_value) != 0 & stats::sd(mean_expression_value) != 0), by = .(locus_name)]
  keep_final_genes <- unique(data_to_transform_df$locus_name[data_to_transform_df$keep_final == TRUE])
  discard_final_genes <- unique(data_to_transform_df$locus_name[data_to_transform_df$keep_final == FALSE])

  mean_df <- mean_df[mean_df$locus_name %in% keep_final_genes, ]
  names(mean_df)[names(mean_df) == "mean_expression_value"] <- "expression_value"

  return(mean_df)
}


#' Get expression of interest
#'
#' @param data_ref File name in working directory, path to file of reference data.
#' @param data_to_transform File name in working directory, path to file of data to transform.
#' @param id_table File name in working directory, path to file of ID table connecting both reference data and data.
#' @param lookup_col_ref_and_id_table Column names shared by both reference data and ID table.
#' @param lookup_col_ref_and_to_transform Column names shared by both reference data and data to transform.
#' @param colnames_wanted List of column names to keep from both reference data and data to transform.
#' @param tissue_wanted Name of tissue from which data will be compared.
#' @param gene_of_interest_acc Gene accession list from data to transform.
#' @param sum_exp_data_ref If \code{TRUE} then sum all gene data. Default is \code{FALSE}.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @return A data frame contains both reference data and data to transform for selected gene of interest.
#' @importFrom rlang .data
#' @export
get_expression_of_interest <- function(data_ref,
                                       data_to_transform,
                                       id_table,
                                       lookup_col_ref_and_id_table = "CDS.model",
                                       lookup_col_ref_and_to_transform = "locus_name",
                                       colnames_wanted = NULL,
                                       tissue_wanted = NULL,
                                       gene_of_interest_acc,
                                       sum_exp_data_ref = FALSE,
                                       accession_data_to_transform = "Col0") {
  # Suppress "no visible binding for global variable" note
  # ..ids_data_ref_colnames <- NULL

  # Load of the single df data
  master_exp <- get_all_data(
    data_ref,
    data_to_transform,
    id_table,
    lookup_col_ref_and_id_table,
    lookup_col_ref_and_to_transform,
    colnames_wanted
  )

  master_exp <- unique(master_exp)

  # Cut down to common tissue if it is specified
  if (is.null(tissue_wanted)) {
    cli::cli_alert_info("No tissue was specified, using all data")
    exp <- master_exp
  } else {
    cli::cli_alert_info("Using data for tissue: {tissue_wanted}")
    exp <- master_exp[master_exp$tissue %in% tissue_wanted, ]
  }

  # Cut down to genes of interest (based on membership of comparison_genes.tsv)
  exp <- exp[exp[[lookup_col_ref_and_to_transform]] %in% gene_of_interest_acc, ]

  # Reformat depending on how want to compare candidate transformed data to fix data, using individual fix data, or summed fix data genes
  # Example: if want to used summed brassica data to compare to the brassica: get symbol level expression total. Locus_name identity is ATG id
  # Otherwise duplicate each arabidopsis, so have an arabidopis copy for each brassica CDS gene. Now locus_name identity is CDS.model
  if (sum_exp_data_ref) {
    exp <- stats::aggregate(
      formula = norm.cpm ~ sample_id + accession + tissue + timepoint + dataset + group + locus_name,
      data = exp,
      FUN = sum
    )
  } else {
    exp_data_to_transform <- exp[exp$accession == accession_data_to_transform, ]
    exp_data_ref <- exp[exp$accession != accession_data_to_transform, ]

    ids_data_ref_colnames <- c(lookup_col_ref_and_id_table, lookup_col_ref_and_to_transform)
    ids_data_ref <- unique(exp_data_ref[, ..ids_data_ref_colnames])
    # exp_data_to_transform$CDS.model <- NULL

    exp_data_to_transform <- exp_data_to_transform %>%
      dplyr::select(-c(lookup_col_ref_and_id_table))

    exp_data_to_transform <- merge(
      ids_data_ref,
      exp_data_to_transform,
      by = lookup_col_ref_and_to_transform,
      allow.cartesian = TRUE
    )

    # Define the same ID and locus name for each data
    exp_data_to_transform <- exp_data_to_transform %>%
      dplyr::mutate(
        id_transform_data = get(lookup_col_ref_and_to_transform),
        locus_name = get(lookup_col_ref_and_id_table)
      )

    exp_data_ref <- exp_data_ref %>%
      dplyr::mutate(
        id_transform_data = get(lookup_col_ref_and_to_transform),
        locus_name = get(lookup_col_ref_and_id_table)
      )

    exp <- rbind(exp_data_to_transform, exp_data_ref)
  }

  return(exp)
}


#' Get all data
#' @noRd
get_all_data <- function(data_ref,
                         data_to_transform,
                         id_table,
                         lookup_col_ref_and_id_table = "CDS.model",
                         lookup_col_ref_and_to_transform = "locus_name",
                         colnames_wanted = NULL) {
  # Suppress "no visible binding for global variable" note
  # ..colnames_id_table <- NULL
  # ..colnames_wanted <- NULL

  # Take unique id_table
  colnames_id_table <- c(lookup_col_ref_and_id_table, lookup_col_ref_and_to_transform)
  id_table_unique <- unique(id_table[, ..colnames_id_table]) %>%
    dplyr::mutate_all(.funs = toupper)

  # Filtering NA values in dataframe of table id
  # dplyr::filter(!is.na(.data$locus_name), !(.data$locus_name %in% c("", "-")))

  # Add fix dataframe info to reg dataframe from
  data_ref <- merge(
    data_ref,
    id_table_unique,
    by = lookup_col_ref_and_id_table
  )

  # Create a column in data_to_transform
  # data_to_transform$locus_name <- data_to_transform$CDS.model
  data_to_transform[, (lookup_col_ref_and_to_transform) := data_to_transform[[lookup_col_ref_and_id_table]]]

  # Cut down to only have genes with ATG locus present in both datasets
  common_symbols <- intersect(
    data_to_transform[[lookup_col_ref_and_to_transform]],
    data_ref[[lookup_col_ref_and_to_transform]]
  )

  data_to_transform <- data_to_transform[data_to_transform[[lookup_col_ref_and_to_transform]] %in% common_symbols, ]
  data_ref <- data_ref[data_ref[[lookup_col_ref_and_to_transform]] %in% common_symbols, ]

  # Take a common columns
  if (is.null(colnames_wanted)) {
    colnames_wanted <- intersect(
      colnames(data_to_transform),
      colnames(data_ref)
    )
  }

  # Join the two datasets into 1 & housekeeping
  expression <- rbind(data_ref[, ..colnames_wanted], data_to_transform[, ..colnames_wanted])

  # Cut down to remove the 'blank' symbol
  expression <- expression[expression[[lookup_col_ref_and_to_transform]] != "", ]

  return(expression)
}
