#' Simplify reference data and data to transform into one dataframe
#'
#' @param filepath_data_ref File name in working directory, path to file of reference data.
#' @param filepath_data_to_transform File name in working directory, path to file of data to transform.
#' @param filepath_id_table File name in working directory, path to file of ID table connecting both reference data and data.
#' @param fix_id_table_shared_colname Column names shared by both reference data and ID table.
#' @param fix_and_to_transform_data_shared_colname Column names shared by both reference data and data to transform.
#' @param colnames_id_table ID table column names.
#' @param colnames_wanted List of column names to keep from both reference data and data to transform.
#' @param tissue_wanted Name of tissue from which data will be compared.
#' @param curr_GoIs Gene of interest list.
#' @param sum_exp_data_ref default is FALSE. If TRUE then sum all gene data.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param ids_data_ref_colnames Column names shared by both reference data and ID table, whose element needs to be unique.
#' @param max_expression_value_wanted Maximum value of expression desired.
#' @param exp_threshold Minimum expression threshold from which expression below the threshold will be removed.
#'
#' @return Combined data frame for both reference data and data to transform: (1) only containing mean of the expression, and (2) still contains replicate data.
get_mean_and_all_exp_data <- function(filepath_data_ref,
                                      filepath_data_to_transform,
                                      filepath_id_table,
                                      fix_id_table_shared_colname,
                                      fix_and_to_transform_data_shared_colname,
                                      colnames_id_table,
                                      colnames_wanted,
                                      tissue_wanted,
                                      curr_GoIs,
                                      sum_exp_data_ref = FALSE,
                                      accession_data_to_transform = "Col0",
                                      ids_data_ref_colnames = c("CDS.model", "locus_name"),
                                      max_expression_value_wanted = 5,
                                      exp_threshold = 0.5) {

  # Load the expression data for all the curr_GoIs gene models, for data to transform and for reference data
  exp <- get_expression_of_interest(
    filepath_data_ref,
    filepath_data_to_transform,
    filepath_id_table,
    fix_id_table_shared_colname,
    fix_and_to_transform_data_shared_colname,
    colnames_id_table,
    colnames_wanted,
    tissue_wanted,
    curr_GoIs,
    sum_exp_data_ref,
    accession_data_to_transform,
    ids_data_ref_colnames
  )

  # Calculate mean of each timepoint by adding a column called "expression_value"
  # TODO: make vector in mean_df a non-hardcoded parameter
  exp[, expression_value := mean(norm.cpm), by = list(locus_name, accession, tissue, timepoint)]
  mean_df <- unique(exp[, c("locus_name", "accession", "tissue", "timepoint", "expression_value")])

  # Filter mean_df to remove genes with very low expression - remove if max is less than 5, and less than half timepoints expressed greater than 1
  data_ref_df <- mean_df[mean_df$accession != accession_data_to_transform]
  data_ref_df[, keep := (max(expression_value) > max_expression_value_wanted | mean(expression_value > 1) > exp_threshold), by = .(locus_name)]

  keep_data_ref_genes <- unique(data_ref_df$locus_name[data_ref_df$keep == TRUE])
  discard_data_ref_genes <- unique(data_ref_df$locus_name[data_ref_df$keep == FALSE])

  # Filter mean_df to remove all data to transform genes with all zeros values
  data_to_transform_df <- mean_df[mean_df$locus_name %in% keep_data_ref_genes & mean_df$accession == accession_data_to_transform]
  data_to_transform_df[, keep_final := (mean(expression_value) != 0 & sd(expression_value) != 0), by = .(locus_name)]
  keep_final_genes <- unique(data_to_transform_df$locus_name[data_to_transform_df$keep_final == TRUE])
  discard_final_genes <- unique(data_to_transform_df$locus_name[data_to_transform_df$keep_final == FALSE])

  mean_df <- mean_df[mean_df$locus_name %in% keep_final_genes, ]

  # Printing the keep genes
  message(paste0(length(keep_data_ref_genes), " brassica genes considered in the comparison"))
  message(paste0(length(keep_final_genes), " all genes considered in the comparison"))

  # Get mean_df, including column "group"
  exp <- exp[exp$locus_name %in% unique(mean_df$locus_name)]
  exp <- subset(exp, select = c("locus_name", "accession", "tissue", "timepoint", "norm.cpm", "group"))
  names(exp)[names(exp) == "norm.cpm"] <- "expression_value"

  # Results object
  results_list <- list(mean_df, exp)

  return(results_list)
}

#' Get expression of interest
get_expression_of_interest <- function(filepath_data_ref,
                                       filepath_data_to_transform,
                                       filepath_id_table,
                                       fix_id_table_shared_colname = "CDS.model",
                                       fix_and_to_transform_data_shared_colname = "locus_name",
                                       colnames_id_table = c("CDS.model", "symbol", "locus_name"),
                                       colnames_wanted = NULL,
                                       tissue_wanted,
                                       curr_GoIs,
                                       sum_exp_data_ref = FALSE,
                                       accession_data_to_transform = "Col0",
                                       ids_data_ref_colnames = c("CDS.model", "locus_name")) {

  # Load of the single df data
  master_exp <- get_all_data(
    filepath_data_ref,
    filepath_data_to_transform,
    filepath_id_table,
    fix_id_table_shared_colname,
    fix_and_to_transform_data_shared_colname,
    colnames_id_table,
    colnames_wanted
  )

  master_exp <- unique(master_exp)

  # Cut down to common tissue if it is specified
  if (is.null(tissue_wanted)) {
    exp <- master_exp
  } else {
    exp <- master_exp[master_exp$tissue %in% tissue_wanted, ]
  }

  # Cut down to genes of interest (based on membership of comparison_genes.tsv)
  exp <- exp[exp$locus_name %in% curr_GoIs, ]

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

    ids_data_ref <- unique(exp_data_ref[, ..ids_data_ref_colnames])
    # exp_data_to_transform$CDS.model <- NULL

    exp_data_to_transform <- exp_data_to_transform %>%
      dplyr::select(-c(fix_id_table_shared_colname))

    exp_data_to_transform <- merge(
      ids_data_ref,
      exp_data_to_transform,
      by = fix_and_to_transform_data_shared_colname,
      allow.cartesian = TRUE
    )

    # Define the same ID and locus name for each data
    exp_data_to_transform <- exp_data_to_transform %>%
      dplyr::mutate(
        id_transform_data = get(fix_and_to_transform_data_shared_colname),
        locus_name = get(fix_id_table_shared_colname)
      )

    exp_data_ref <- exp_data_ref %>%
      dplyr::mutate(
        id_transform_data = get(fix_and_to_transform_data_shared_colname),
        locus_name = get(fix_id_table_shared_colname)
      )

    # exp_data_to_transform$ara.id <- exp_data_to_transform$locus_name
    # exp_data_to_transform$locus_name <- exp_data_to_transform$CDS.model
    # exp_data_ref$ara.id <- exp_data_ref$locus_name
    # exp_data_ref$locus_name <- exp_data_ref$CDS.model

    exp <- rbind(exp_data_to_transform, exp_data_ref)
  }

  # Shorten experiment group names
  exp <- shorten_groups(exp)
  return(exp)
}


#' Get all data
get_all_data <- function(filepath_data_ref,
                         filepath_data_to_transform,
                         filepath_id_table,
                         fix_id_table_shared_colname = "CDS.model",
                         fix_and_to_transform_data_shared_colname = "locus_name",
                         colnames_id_table = c("CDS.model", "symbol", "locus_name"),
                         colnames_wanted = NULL) {

  # Read RDS file
  data_ref <- readRDS(filepath_data_ref)
  data_to_transform <- readRDS(filepath_data_to_transform)

  if (tools::file_ext(filepath_id_table) == "csv") {
    id_table <- data.table::fread(filepath_id_table)
  } else {
    id_table <- readRDS(filepath_id_table)
  }

  # Take unique id_table
  id_table_unique <- unique(id_table[, ..colnames_id_table]) %>%
    dplyr::mutate_all(.funs = toupper)

  # Filtering NA values in dataframe of table id
  # dplyr::filter(!is.na(locus_name), !locus_name %in% c("", "-"))

  # Add fix dataframe info to reg dataframe from
  data_ref <- merge(data_ref, id_table_unique, by = fix_id_table_shared_colname)

  # Create a column in data_to_transform
  # data_to_transform$locus_name <- data_to_transform$CDS.model
  data_to_transform[, (fix_and_to_transform_data_shared_colname) := data_to_transform[[fix_id_table_shared_colname]]]

  # Cut down to only have genes with ATG locus present in both datasets
  common_symbols <- intersect(
    data_to_transform[[fix_and_to_transform_data_shared_colname]],
    data_ref[[fix_and_to_transform_data_shared_colname]]
  )

  data_to_transform <- data_to_transform[data_to_transform[[fix_and_to_transform_data_shared_colname]] %in% common_symbols, ]
  data_ref <- data_ref[data_ref[[fix_and_to_transform_data_shared_colname]] %in% common_symbols, ]

  # Take a common columns
  if (is.null(colnames_wanted)) {
    colnames_wanted <- intersect(
      colnames(data_to_transform),
      colnames(data_ref)
    )
  } else {
    colnames_wanted <- colnames_wanted
  }

  # Join the two datasets into 1 & housekeeping
  expression <- rbind(data_ref[, ..colnames_wanted], data_to_transform[, ..colnames_wanted])

  # Cut down to remove the 'blank' symbol
  expression <- expression[expression[[fix_and_to_transform_data_shared_colname]] != "", ]

  return(expression)
}


#' Shorten groups
#' @param exp
#' @param accession_data_to_transform
shorten_groups <- function(exp) {
  # Get reps for klepikova and for brassica data, make sure it is a data.table
  exp <- data.table::data.table(exp)

  # Get the last element of "sample_id" of Brassica data
  B <- exp[exp$accession != "Col0"]
  B[, c("j1", "j2", "j3", "j4", "j5", "j6", "rep") := data.table::tstrsplit(sample_id, split = "_")]
  B$rep[is.na(B$rep)] <- 1
  B[, c("j1", "j2", "j3", "j4", "j5", "j6")] <- NULL

  # Get the last element of "dataset" of Arabidopsis data
  A <- exp[exp$accession == "Col0"]
  A[, c("j1", "j2", "rep") := data.table::tstrsplit(dataset, split = "_")]
  A[, c("j1", "j2")] <- NULL
m
  # Combine both data
  exp <- rbind(B, A)

  exp$ds <- factor(
    exp$rep,
    levels = c("1", "2", "3", "4"),
    labels = c("a", "b", "c", "d")
  ) %>%
    as.character()

  exp <- data.table::data.table(exp)
  exp[, group := paste(accession, sprintf("%02d", timepoint), ds, sep = "-")]

  return(exp)
}
