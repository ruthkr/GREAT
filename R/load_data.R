#' @export
load_mean_df <- function(filepath_data_target,
                          filepath_data_to_transform,
                          filepath_id_table,
                          target_id_table_shared_colname,
                          target_and_to_transform_data_shared_colname,
                          colnames_id_table,
                          colnames_wanted,
                          tissue_wanted,
                          curr_GoIs,
                          sum_exp_data_target,
                          accession_data_to_transform,
                          ids_data_target_colnames,
                          max_mean_cpm_wanted = 5) {

  # Load the expression data for all the curr_GoIs gene models, for data to transform and for data target
  exp <- get_expression_of_interest(filepath_data_target,
                                    filepath_data_to_transform,
                                    filepath_id_table,
                                    target_id_table_shared_colname,
                                    target_and_to_transform_data_shared_colname,
                                    colnames_id_table,
                                    colnames_wanted,
                                    tissue_wanted,
                                    curr_GoIs,
                                    sum_exp_data_target = F,
                                    accession_data_to_transform = "Col0",
                                    ids_data_target_colnames = c("CDS.model", "locus_name"))

  # Calculate mean of each timepoint by adding a column called "mean_cpm"
  exp[, mean_cpm := mean(norm.cpm), by = list(locus_name, accession, tissue, timepoint)]
  mean_df <- unique(exp[, c('locus_name', 'accession', 'tissue', 'timepoint', 'mean_cpm')])

  # Filter mean_df to remove genes with very low expression - remove if max is less than 5, and less than half timepoints expressed greater than 1
  data_target_df <- mean_df[mean_df$accession != accession_data_to_transform]
  data_target_df[, keep := (max(mean_cpm) > max_mean_cpm_wanted | mean(mean_cpm > 1) > 0.5) , by = .(locus_name)]

  keep_data_target_genes <- unique(data_target_df$locus_name[data_target_df$keep == TRUE])
  discard_data_target_genes <- unique(data_target_df$locus_name[data_target_df$keep == FALSE])

  # Filter mean_df to remove all data to transform genes with all zeros values
  data_to_transform_df <- mean_df[mean_df$locus_name %in% keep_data_target_genes & mean_df$accession == accession_data_to_transform]
  data_to_transform_df[, keep_final:=(mean(mean_cpm) != 0 & sd(mean_cpm) != 0), by=.(locus_name)]
  keep_final_genes <- unique(data_to_transform_df$locus_name[data_to_transform_df$keep_final==TRUE])
  discard_final_genes <- unique(data_to_transform_df$locus_name[data_to_transform_df$keep_final==FALSE])

  mean_df <- mean_df[mean_df$locus_name %in% keep_final_genes, ]

  # Printing the keep genes
  print(paste0(length(keep_data_target_genes), ' brassica genes considered in the comparison'))
  print(paste0(length(keep_final_genes), ' all genes considered in the comparison'))

  # Get mean_df, including column "group"
  exp <- exp[exp$locus_name %in% unique(mean_df$locus_name)]
  exp <- subset(exp, select = c('locus_name', 'accession', 'tissue', 'timepoint',
                              'norm.cpm', 'group'))
  names(exp)[names(exp) == 'norm.cpm'] <- 'mean_cpm'
  return(list(mean_df, exp))
}

#' @export
get_expression_of_interest <- function(filepath_data_target,
                                       filepath_data_to_transform,
                                       filepath_id_table,
                                       target_id_table_shared_colname,
                                       target_and_to_transform_data_shared_colname,
                                       colnames_id_table,
                                       colnames_wanted,
                                       tissue_wanted,
                                       curr_GoIs,
                                       sum_exp_data_target = F,
                                       accession_data_to_transform = "Col0",
                                       ids_data_target_colnames = c("CDS.model", "locus_name")) {

  # Load of the single df data
  master_exp <- get_all_data(filepath_data_target,
                           filepath_data_to_transform,
                           filepath_id_table,
                           target_id_table_shared_colname = "CDS.model",
                           target_and_to_transform_data_shared_colname = "locus_name",
                           colnames_id_table = c("CDS.model", "symbol", "locus_name"),
                           colnames_wanted = NULL)

  master_exp <- unique(master_exp)

  # Cut down to common tissue
  exp <- master_exp[master_exp$tissue %in% tissue_wanted, ]

  # Cut down to genes of interest (based on membership of comparison_genes.tsv)
  exp <- exp[exp$locus_name %in% curr_GoIs, ]

  # Reformat depending on how want to compare candidate transformed data to target data, using individual target data, or summed target data genes
  # Example: if want to used summed brassica data to compare to the brassica: get symbol level expression total. Locus_name identity is ATG id
  # Otherwise duplicate each arabidopsis, so have an arabidopis copy for each brassica CDS gene. Now locus_name identity is CDS.model
  if (sum_exp_data_target == T) {
    exp <- stats::aggregate(norm.cpm~sample_id+accession+tissue+timepoint+dataset+group+locus_name, data=exp, sum)
  } else if (sum_exp_data_target == F) {
    exp_data_to_transform <- exp[exp$accession == accession_data_to_transform,]
    exp_data_target <- exp[exp$accession != accession_data_to_transform,]

    ids_data_target <- unique(exp_data_target[, ..ids_data_target_colnames])
    # exp_data_to_transform$CDS.model <- NULL

    exp_data_to_transform <- exp_data_to_transform %>%
      dplyr::select(-c(target_id_table_shared_colname))

    exp_data_to_transform <- merge(ids_data_target, exp_data_to_transform, by = target_and_to_transform_data_shared_colname, allow.cartesian=T)

    # Define the same ID and locus name for each data
    exp_data_to_transform <- exp_data_to_transform %>%
      dplyr::mutate(
        id_transform_data = get(target_and_to_transform_data_shared_colname),
        locus_name = get(target_id_table_shared_colname)
      )

    exp_data_target <- exp_data_target %>%
      dplyr::mutate(
        id_transform_data = get(target_and_to_transform_data_shared_colname),
        locus_name = get(target_id_table_shared_colname)
      )


    # exp_data_to_transform$ara.id <- exp_data_to_transform$locus_name
    # exp_data_to_transform$locus_name <- exp_data_to_transform$CDS.model
    # exp_data_target$ara.id <- exp_data_target$locus_name
    # exp_data_target$locus_name <- exp_data_target$CDS.model


    exp <- rbind(exp_data_to_transform, exp_data_target)
  }

  # Shorten experiment group names
  exp <- shorten_groups(exp)
  return(exp)
}


#' @export
get_all_data <- function(filepath_data_target,
                         filepath_data_to_transform,
                         filepath_id_table,
                         target_id_table_shared_colname = "CDS.model",
                         target_and_to_transform_data_shared_colname = "locus_name",
                         colnames_id_table = c("CDS.model", "symbol", "locus_name"),
                         colnames_wanted = NULL) {

  # Read RDS file
  data_target <- readRDS(filepath_data_target)
  data_to_transform <- readRDS(filepath_data_to_transform)

  if (tools::file_ext(filepath_id_table) == "csv"){
    id_table <- data.table::fread(filepath_id_table)
  } else {
    id_table <- readRDS(filepath_id_table)
  }

  # Take unique id_table
  id_table_unique <- unique(id_table[, ..colnames_id_table]) %>%
    dplyr:: mutate_all(.funs = toupper)

  # Filtering NA values in dataframe of table id
  # dplyr::filter(!is.na(locus_name), !locus_name %in% c("", "-"))

  # Add target dataframe info to reg dataframe from
  data_target <- merge(data_target, id_table_unique, by = target_id_table_shared_colname)

  # Create a column in data_to_transform
  # data_to_transform$locus_name <- data_to_transform$CDS.model
  data_to_transform[, (target_and_to_transform_data_shared_colname) := data_to_transform[[target_id_table_shared_colname]]]

  # Cut down to only have genes with ATG locus present in both datasets
  common_symbols <- intersect(data_to_transform[[target_and_to_transform_data_shared_colname]], data_target[[target_and_to_transform_data_shared_colname]])
  data_to_transform <- data_to_transform[data_to_transform[[target_and_to_transform_data_shared_colname]] %in% common_symbols, ]
  data_target <- data_target[data_target[[target_and_to_transform_data_shared_colname]] %in% common_symbols, ]

  # Take a common columns
  if (is.null(colnames_wanted)) {
    colnames_wanted <- intersect(colnames(data_to_transform), colnames(data_target))
  } else {
    colnames_wanted <- colnames_wanted
  }

  # Join the two datasets into 1 & housekeeping
  expression <- rbind(data_target[, ..colnames_wanted], data_to_transform[, ..colnames_wanted])

  # Cut down to remove the 'blank' symbol
  expression <- expression[expression[[target_and_to_transform_data_shared_colname]] != "", ]

  return(expression)
}

#' @export
shorten_groups <- function(exp) {

  # Get reps for klepikova and for brassica data, make sure it is a data.table
  exp <- data.table::data.table(exp)

  # Get the last element of "sample_id" of Brassica data
  B <- exp[exp$accession != 'Col0']
  B[, c('j1', 'j2', 'j3', 'j4', 'j5', 'j6', 'rep'):=data.table::tstrsplit(sample_id, split='_')]
  B$rep[is.na(B$rep)] <- 1
  B[, c('j1', 'j2', 'j3', 'j4', 'j5', 'j6')] <- NULL

  # Get the last element of "dataset" of Arabidopsis data
  A <- exp[exp$accession == 'Col0']
  A[, c('j1', 'j2','rep'):=data.table::tstrsplit(dataset, split='_')]
  A[, c('j1', 'j2')] <- NULL

  # Combine both data
  exp <- rbind(B,A)

  exp$ds <- as.character(factor(exp$rep, levels=c('1', '2', '3', '4'), labels=c('a', 'b', 'c', 'd')))

  exp <- data.table::data.table(exp)
  exp[, group:=paste(accession, sprintf('%02d', timepoint), ds, sep='-')]
  return(exp)
}
