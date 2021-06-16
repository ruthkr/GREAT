#' @export
load_mean_df_copy <- function(file_path_brassica, file_path_arabidopsis, file_path_id_table, tissue_wanted, curr_GoIs, sum_brassicas = F) {

  # Load the expression data for all the curr_GoIs gene models, for arabidopsis, and for the specified brassica
  exp <- get_expression_of_interest_copy(file_path_brassica, file_path_arabidopsis, file_path_id_table, tissue_wanted, curr_GoIs, sum_brassicas = F)

  # Calculate mean of each timepoint by adding a column called "mean.cpm"
  exp[, mean.cpm:=mean(norm.cpm), by=list(locus_name, accession, tissue, timepoint)]
  mean_df <- unique(exp[, c('locus_name', 'accession', 'tissue', 'timepoint', 'mean.cpm')])

  # Filter mean_df to remove genes with very low expression - remove if max is less than 5, and less than half timepoints expressed greater than 1
  # bra_df <- mean_df[mean_df$accession != 'Col0']
  # bra_df[, keep:=(max(mean.cpm) > 5 | mean(mean.cpm > 1) > 0.5) , by=.(locus_name)]
  # keep.genes <- unique(bra_df$locus_name[bra_df$keep==TRUE])
  # discard.genes <- unique(bra_df$locus_name[bra_df$keep==FALSE])
  bra_df <- mean_df[mean_df$accession != 'Col0']
  bra_df[, keep:=(max(mean.cpm) > 5 | mean(mean.cpm > 1) > 0.5) , by=.(locus_name)]
  keep_bra_genes <- unique(bra_df$locus_name[bra_df$keep==TRUE])
  discard_bra_genes <- unique(bra_df$locus_name[bra_df$keep==FALSE])

  # Filter mean_df to remove all arabidopsis genes with all zeros values
  ara_df <- mean_df[mean_df$locus_name %in% keep_bra_genes & mean_df$accession == 'Col0']
  ara_df[, keep_final:=(mean(mean.cpm) != 0 & sd(mean.cpm) != 0), by=.(locus_name)]
  keep_final_genes <- unique(ara_df$locus_name[ara_df$keep_final==TRUE])
  discard_final_genes <- unique(ara_df$locus_name[ara_df$keep_final==FALSE])

  mean_df <- mean_df[mean_df$locus_name %in% keep_final_genes,]

  # Printing the keep genes
  print(paste0(length(keep_bra_genes), ' brassica genes considered in the comparison'))
  print(paste0(length(keep_final_genes), ' all genes considered in the comparison'))


  # print(paste0(length(unique(mean_df$locus_name)), ' brassica genes considered in the comparison'))
  # print(paste(c("Discarded genes:", paste(discard.genes, collapse = ", ")), collapse = " "))

  # Get mean_df, including column "group"
  exp <- exp[exp$locus_name %in% unique(mean_df$locus_name)]
  exp <- subset(exp, select=c('locus_name', 'accession', 'tissue', 'timepoint',
                              'norm.cpm', 'group'))
  names(exp)[names(exp)=='norm.cpm'] <- 'mean.cpm'
  return(list(mean_df, exp))
}

#' @export
get_expression_of_interest_copy <- function(file_path_brassica, file_path_arabidopsis, file_path_id_table, tissue_wanted, curr_GoIs, sum_brassicas = F) {

  # Load rds and arabidopsis gene expression data into single df.
  master_exp <- get_all_data_copy(file_path_brassica, file_path_arabidopsis, file_path_id_table)
  master_exp <- unique(master_exp)

  # Cut down to common tissue
  exp <- master_exp[master_exp$tissue == tissue_wanted, ]

  # Cut down to genes of interest (based on membership of comparison_genes.tsv)
  exp <- exp[exp$locus_name %in% curr_GoIs,]

  # Reformat depending on how want to compare arabidopsis to brassica, using indiv. brassica genes, or summed brassica genes
  # if want to used summed brassica data to compare to the brassica: get symbol level expression total. Locus_name identity is ATG id
  if (sum_brassicas == T) {
    exp <- stats::aggregate(norm.cpm~sample_id+accession+tissue+timepoint+dataset+group+locus_name, data=exp, sum)
  } else if (sum_brassicas == F) {
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
get_all_data_copy <- function(file_path_brassica, file_path_arabidopsis, file_path_id_table, colnames_wanted = NULL) {

  # Read RDS file
  bra_data <- readRDS(file_path_brassica)
  ara_data <- readRDS(file_path_arabidopsis)

  if (tools::file_ext(file_path_id_table) == "csv"){
    id_table <- data.table::fread(file_path_id_table)
  } else {
    id_table <- readRDS(file_path_id_table)
  }


  id_table <- unique(id_table[, c("locus_name", "symbol", "CDS.model")])

  # Take unique id_table
  id_table_unique <- unique(id_table[, c("CDS.model", "symbol", "locus_name")]) %>%
    dplyr::mutate(CDS.model = toupper(CDS.model)) %>%
    dplyr::filter(!is.na(locus_name), !locus_name %in% c("", "-"))

  # Add ATG locus info
  bra_data <- merge(bra_data, id_table_unique, by = "CDS.model")

  # Create a column in ara_data
  ara_data$locus_name <- ara_data$CDS.model

  # cut down to only have genes with ATG locus present in both datasets
  common_symbols <- intersect(ara_data$locus_name, bra_data$locus_name)
  ara_data <- ara_data[ara_data$locus_name %in% common_symbols, ]
  bra_data <- bra_data[bra_data$locus_name %in% common_symbols, ]

  # Take a common columns
  if (is.null(colnames_wanted)) {
    colnames_wanted <- intersect(colnames(ara_data), colnames(bra_data))
  } else {
    colnames_wanted <- colnames_wanted
  }


  # Join the two datasets into 1 & housekeeping
  expression <- rbind(bra_data[, ..colnames_wanted], ara_data[, ..colnames_wanted])

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
  B[, c('j1', 'j2', 'j3', 'j4', 'j5', 'j6', 'rep') := data.table::tstrsplit(sample_id, split='_')]
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
