#' @export
load_mean_df <- function(filepath_reg_data, filepath_target_data, filepath_id_table, tissue_wanted, curr_GoIs, sum_exp_reg_data = F) {

  # Load the expression data for all the curr_GoIs gene models, for arabidopsis, and for the specified brassica
  exp <- get_expression_of_interest(filepath_reg_data, filepath_target_data, filepath_id_table, tissue_wanted, curr_GoIs, sum_exp_reg_data = F)

  # Calculate mean of each timepoint by adding a column called "mean.cpm"
  exp[, mean.cpm:=mean(norm.cpm), by=list(locus_name, accession, tissue, timepoint)]
  mean.df <- unique(exp[, c('locus_name', 'accession', 'tissue', 'timepoint', 'mean.cpm')])

  # Filter mean.df to remove genes with very low expression - remove if max is less than 5, and less than half timepoints expressed greater than 1
  # bra_df <- mean.df[mean.df$accession != 'Col0']
  # bra_df[, keep:=(max(mean.cpm) > 5 | mean(mean.cpm > 1) > 0.5) , by=.(locus_name)]
  # keep.genes <- unique(bra_df$locus_name[bra_df$keep==TRUE])
  # discard.genes <- unique(bra_df$locus_name[bra_df$keep==FALSE])
  bra_df <- mean.df[mean.df$accession != 'Col0']
  bra_df[, keep:=(max(mean.cpm) > 5 | mean(mean.cpm > 1) > 0.5) , by=.(locus_name)]
  keep_bra_genes <- unique(bra_df$locus_name[bra_df$keep==TRUE])
  discard_bra_genes <- unique(bra_df$locus_name[bra_df$keep==FALSE])

  # Filter mean.df to remove all arabidopsis genes with all zeros values
  ara_df <- mean.df[mean.df$locus_name %in% keep_bra_genes & mean.df$accession == 'Col0']
  ara_df[, keep_final:=(mean(mean.cpm) != 0 & sd(mean.cpm) != 0), by=.(locus_name)]
  keep_final_genes <- unique(ara_df$locus_name[ara_df$keep_final==TRUE])
  discard_final_genes <- unique(ara_df$locus_name[ara_df$keep_final==FALSE])

  mean.df <- mean.df[mean.df$locus_name %in% keep_final_genes,]

  # Printing the keep genes
  print(paste0(length(keep_bra_genes), ' brassica genes considered in the comparison'))
  print(paste0(length(keep_final_genes), ' all genes considered in the comparison'))


  # print(paste0(length(unique(mean.df$locus_name)), ' brassica genes considered in the comparison'))
  # print(paste(c("Discarded genes:", paste(discard.genes, collapse = ", ")), collapse = " "))

  # Get mean.df, including column "group"
  exp <- exp[exp$locus_name %in% unique(mean.df$locus_name)]
  exp <- subset(exp, select=c('locus_name', 'accession', 'tissue', 'timepoint',
                              'norm.cpm', 'group'))
  names(exp)[names(exp)=='norm.cpm'] <- 'mean.cpm'
  return(list(mean.df, exp))
}

#' @export
get_expression_of_interest <- function(filepath_reg_data, filepath_target_data, filepath_id_table, tissue_wanted, curr_GoIs, sum_exp_reg_data = F) {

  # Load rds and arabidopsis gene expression data into single df.
  master_exp <- get_all_data(filepath_reg_data, filepath_target_data, filepath_id_table)
  master_exp <- unique(master_exp)

  # Cut down to common tissue
  exp <- master_exp[master_exp$tissue == tissue_wanted, ]

  # Cut down to genes of interest (based on membership of comparison_genes.tsv)
  exp <- exp[exp$locus_name %in% curr_GoIs,]

  # Reformat depending on how want to compare arabidopsis to brassica, using indiv. brassica genes, or summed brassica genes
  # if want to used summed brassica data to compare to the brassica: get symbol level expression total. Locus_name identity is ATG id
  if (sum_exp_reg_data == T) {
    exp <- stats::aggregate(norm.cpm~sample_id+accession+tissue+timepoint+dataset+group+locus_name, data=exp, sum)
  } else if (sum_exp_reg_data == F) {
    # Otherwise duplicate each arabidopsis, so have an arabidopis copy for each brassica CDS gene. Now locus_name identity is CDS.model
    ara.exp <- exp[exp$accession=='Col0',]
    bra.exp <- exp[exp$accession!='Col0',]
    bra.ids <- unique(bra.exp[, c('CDS.model', 'locus_name')])
    ara.exp$CDS.model <- NULL
    ara.exp <- merge(bra.ids, ara.exp, by='locus_name', allow.cartesian=T)
    ara.exp$ara.id <- ara.exp$locus_name
    ara.exp$locus_name <- ara.exp$CDS.model
    bra.exp$ara.id <- bra.exp$locus_name
    bra.exp$locus_name <- bra.exp$CDS.model
    exp <- rbind(ara.exp, bra.exp)
  }

  # Shorten experiment group names
  exp <- shorten_groups(exp)
  return(exp)
}


#' @export
get_all_data <- function(filepath_reg_data, filepath_target_data, filepath_id_table, colnames_id_table = c("CDS.model", "symbol", "locus_name"), colnames_wanted = NULL) {

  # Read RDS file
  reg_data <- readRDS(filepath_reg_data)
  target_data <- readRDS(filepath_target_data)

  if (tools::file_ext(filepath_id_table) == "csv"){
    id_table <- data.table::fread(filepath_id_table)
  } else {
    id_table <- readRDS(filepath_id_table)
  }

  # Take unique id_table
  id_table_unique <- unique(id_table[, colnames_id_table]) %>%
    dplyr::mutate(CDS.model = toupper(CDS.model)) %>%
    dplyr::filter(!is.na(locus_name), !locus_name %in% c("", "-"))

  # Add ATG locus info
  reg_data <- merge(reg_data, id_table_unique, by = "CDS.model")

  # Create a column in target_data
  target_data$locus_name <- target_data$CDS.model

  # cut down to only have genes with ATG locus present in both datasets
  common_symbols <- intersect(target_data$locus_name, reg_data$locus_name)
  target_data <- target_data[target_data$locus_name %in% common_symbols, ]
  reg_data <- reg_data[reg_data$locus_name %in% common_symbols, ]

  # Take a common columns
  if (is.null(colnames_wanted)) {
    colnames_wanted <- intersect(colnames(target_data), colnames(reg_data))
  } else {
    colnames_wanted <- colnames_wanted
  }


  # Join the two datasets into 1 & housekeeping
  expression <- rbind(reg_data[, ..colnames_wanted], target_data[, ..colnames_wanted])

  # Cut down to remove the 'blank' symbol
  expression <- expression[expression$locus_name != "", ]

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
