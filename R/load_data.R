#' @export
load_mean.df <- function() {

  #setwd('/Volumes/Research-Projects/bravo/alex/BRAVO_rna-seq/scripts/')
  rds_file <- 'ro18_chiifu_apex' # don't include the .rds # the name of the brassica data to load
  sumBrassicas <- F # if false use seperate brassica genes, and compare to repeated arabidopsis genes. If true, sume copies of each brassica and compare to arabidopsis
  #datapath <- paste0('../final_data/rds/', rds_file, '.rds')


  #### specify the genes to be used in the comparison: ---
  # GoIs <- data.table::fread(paste0('../graphs/', rds_file, '/comparison_genes.tsv')) # read the list of Arabidopsis id's for the genes of interest

  # ------ changed by Ruth
  GoIs <- data.table::fread(paste0('graphs/', rds_file, '/comparison_genes.tsv'))

  print(unique(GoIs$model))
  keep_model_set <- c('rf') #, 'logistic', 'trough', 'spike', 'linear') # rf means don't force them to be the same AT ALL!
  curr_GoIs <- GoIs$symbol[GoIs$model %in% keep_model_set]
  # add some particular flowering genes of interest
  # AT2G22540 : SVP
  # AT3G57920 : SPL15
  # AT2G42200 : SPL9
  # AT5G51870 : AGL71 : BRAA10G010470.3C
  # AT5G51860 : AGL72
  # AT4G36920 : AP2
  curr_GoIs <- c(curr_GoIs, 'AT2G22540', 'AT3G57920', 'AT2G42200', 'AT5G51870', 'AT5G51860', 'AT4G36920')
  #curr_GoIs <- c('LMI1', 'STM', 'PNY', 'PNF', 'UFO', 'AGL42', 'AGL71', 'AGL72', 'PI', 'GA20OX1', 'GA20OX2') # KEY THINGS, MISSING



  ##### load the expression data for all the curr_GoIs gene models, for arabidopsis, and for the specified brassica
  filt_models <- get_expression_oI(rds_file, curr_GoIs, sumBrassicas) #
  length(unique(curr_GoIs))
  length(unique(filt_models$CDS.model))
  exp <- filt_models

  # get mean of each timepoint
  exp[, mean.cpm:=mean(norm.cpm), by=list(locus_name, accession, tissue, timepoint)]
  mean.df <- unique(exp[, c('locus_name', 'accession', 'tissue', 'timepoint', 'mean.cpm')])



  # filter mean.df to remove genes with very low expression - remove if max is less than 5, and less than half timepoints expressed
  # greater than 1
  ro18.df <- mean.df[mean.df$accession=='Ro18']
  ro18.df[, keep:=(max(mean.cpm) > 2 | mean(mean.cpm > 1) > 0.5) , by=.(locus_name)]
  keep.genes <- unique(ro18.df$locus_name[ro18.df$keep==TRUE])
  mean.df <- mean.df[mean.df$locus_name %in% keep.genes,]
  print(paste0(length(unique(mean.df$locus_name)), ' genes considered in the comparison'))
  rm(ro18.df, keep.genes)

  exp <- exp[exp$locus_name %in% unique(mean.df$locus_name)]
  exp <- subset(exp, select=c('locus_name', 'accession', 'tissue', 'timepoint',
                              'norm.cpm', 'group'))
  names(exp)[names(exp)=='norm.cpm'] <- 'mean.cpm'
  return(list(mean.df, exp))
}

#' @export
get_expression_oI <- function(rds_file, curr_GoIs, sumBrassicas) {

  # load rds and arabidopsis gene expression data into single df.
  master_exp <- get_all_data(rds_file)
  master_exp <- unique(master_exp)
  # cut down to common tissue
  exp <- master_exp[master_exp$tissue=='apex', ]
  #cut down to genes of interest (based on membership of comparison_genes.tsv)
  exp <- exp[exp$locus_name %in% curr_GoIs,]
  tmp <- unique(exp)
  tmp <- exp[exp$CDS.model=='BRAA06G036760.3C',]

  # different models can potentially have selected the same genes, but with diff symbol names - e.g. AP1, ATAP1 - want to ensure only have one of each
  # of these.
  # exp <- exp[order(exp$symbol, decreasing=F), ]
  # # get rid of different symbols for same genes really!
  # tmp <- exp[accession=='Col0', c('CDS.model', 'symbol')]
  # unique_tmp <- tmp[!duplicated(tmp[, c('CDS.model')])] # when multiple symbols for same CDS, keep one arbitrary one.
  # exp <- exp[exp$symbol %in% unique_tmp$symbol,]


  #length(unique(exp$symbol))
  #exp <- data.frame(exp)
  #exp <- exp[!duplicated(exp[c('sample_id', 'CDS.model')]), ]
  #length(unique(test$symbol))

  # reformat depending on how want to compare arabidopsis to brassica, using indiv. brassica genes, or summed brassica genes
  # if want to used summed brassica data to compare to the brassica: get symbol level expression total. Locus_name identity is ATG id
  if (sumBrassicas==T) {
    exp <- stats::aggregate(norm.cpm~sample_id+accession+tissue+timepoint+dataset+group+locus_name, data=exp, sum)
  } else if (sumBrassicas==F) {
    # otherwise duplicate each arabidopsis, so have an arabidopis copy for each brassica CDS gene. Now locus_name identity is CDS.model
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

  # shorten experiment group names
  exp <- shorten_groups(exp)
  return(exp)
}


#' @export
get_all_data <- function(file_path_brassica, file_path_arabidopsis, file_path_id_table) {

  # Read RDS file
  bra_data <- readRDS(file_path_brassica)
  ara_data <- readRDS(file_path_arabidopsis)
  id_table <- readRDS(file_path_id_table)

  id_table <- unique(id_table[, c("locus_name", "symbol", "CDS.model")])

  # Take unique id_table
  id_table_unique <- unique(id_table[, c("CDS.model", "locus_name")]) %>%
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
  colnames_intersect <- intersect(colnames(ara_data), colnames(bra_data))

  # Join the two datasets into 1 & housekeeping
  expression <- rbind(bra_data[, ..colnames_intersect], ara_data[, ..colnames_intersect])

  # Cut down to remove the 'blank' symbol
  expression <- expression[expression$locus_name != "", ]

  return(expression)
}

#' @export
shorten_groups <- function(exp) {
  # get reps for klepikova and for brassica data
  exp <- data.table::data.table(exp)
  B <- exp[exp$accession != 'Col0']
  B[, c('j1', 'j2', 'j3', 'j4', 'j5', 'j6', 'rep'):=data.table::tstrsplit(sample_id, split='_')]
  B$rep[is.na(B$rep)] <- 1
  B[, c('j1', 'j2', 'j3', 'j4', 'j5', 'j6')] <- NULL

  A <- exp[exp$accession == 'Col0']
  A[, c('j1', 'j2','rep'):=data.table::tstrsplit(dataset, split='_')]
  A[, c('j1', 'j2')] <- NULL

  exp <- rbind(B,A)

  # exp$ds <- plyr::mapvalues(exp$rep, from=c('1', '2', '3', '4'), to=c('a', 'b', 'c', 'd'))
  exp$ds <- as.character(factor(exp$rep, levels=c('1', '2', '3', '4'), labels=c('a', 'b', 'c', 'd')))

  exp <- data.table::data.table(exp)
  exp[, group:=paste(accession, sprintf('%02d', timepoint), ds, sep='-')]
  return(exp)
}
