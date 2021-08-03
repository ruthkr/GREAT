# mean_df <- real.mean_df
# mean_df.sc <- real.sc.df
# imputed.mean_df <- imputed.mean_df
#' @export
calculate_between_sample_distance <- function(mean_df, mean_df.sc, imputed.mean_df) {

  ### convert all to wide format ready for distance calculation

  # mean_df
  sample.id.cols <- c('accession','timepoint')
  gene.col <- c('locus_name')
  expression.col <- 'mean_cpm'
  mean.dt.w <- reformat_for_distance_calculation(mean_df, sample.id.cols, gene.col, expression.col)

  # normalised mean_df
  sample.id.cols <- c('accession','timepoint')
  gene.col <- c('locus_name')
  expression.col <- c('sc.mean_cpm')
  mean.dt.sc.w <- reformat_for_distance_calculation(mean_df.sc, sample.id.cols, gene.col, expression.col)

  # imputed.mean_df - all genes
  sample.id.cols <- c('accession','shifted_time')
  gene.col <- c('locus_name')
  expression.col <- c('mean_cpm')
  imputed.mean.dt.w <- reformat_for_distance_calculation(imputed.mean_df, sample.id.cols, gene.col, expression.col)

  # same, but for subsets of REGISTERED / NOT REGISTERED genes.
  # distance between samples, only using genes which are found best model is not registered
  not.registered.genes <- unique(imputed.mean_df$locus_name[imputed.mean_df$is.registered==FALSE])
  mean.dt.sc.w.not.registered <- mean.dt.sc.w[mean.dt.sc.w$locus_name %in% not.registered.genes,]
  # distance between samples, only using genes which are found best when ARE registered
  registered.genes <- unique(imputed.mean_df$locus_name[imputed.mean_df$is.registered==TRUE])
  mean.dt.sc.w.registered <- mean.dt.sc.w[mean.dt.sc.w$locus_name %in% registered.genes,]
  # after registration, but only for registered genes
  imputed.mean.dt.w.registered <- imputed.mean.dt.w[imputed.mean.dt.w$locus_name %in% registered.genes,]


  ### calculate distance between each sample in each of the wide format tables.
  # shifting genes results in NAs for some genes at some timepoints.
  # therefore different numbers of dimensions (genes) depending on the comparison.
  # therefore euclidean distance is not appropriate.
  # Distance used is mean of [squared distance between each gene / absolute mean of expression
  # of that gene in the sample]. Is calculated for each gene for which have data in both samples.
  # Mean of these values (divided by number of genes is calculated for) is reported.
  D.mean <- calc_sample_distance(mean.dt.w, gene.col='locus_name')
  D.scaled <- calc_sample_distance(mean.dt.sc.w, gene.col='locus_name')
  D.registered <- calc_sample_distance(imputed.mean.dt.w, gene.col='locus_name')

  D.scaled.not.registered.genes <- calc_sample_distance(mean.dt.sc.w.not.registered, gene.col='locus_name')
  D.scaled.registered.genes <- calc_sample_distance(mean.dt.sc.w.registered, gene.col='locus_name')
  D.registered.registered.genes <- calc_sample_distance(imputed.mean.dt.w.registered, gene.col='locus_name')


  # for use to make heatmaps with shared scales
  D.mean$title <- 'mean expression'
  D.scaled$title <- 'scaled mean expression (all genes)'
  D.registered$title <- 'registered & scaled mean expression (all genes)'

  D.scaled.not.registered.genes$title <- 'scaled mean expression (only not-registered genes)'
  D.scaled.registered.genes$title <- 'scaled mean expression (only registered genes)'
  D.registered.registered.genes$title <- 'registered & scaled mean expression (only registered genes)'


  return(list('D.mean'=D.mean,
              'D.scaled'=D.scaled,
              'D.registered'=D.registered,
              'D.scaled.onlyNR'=D.scaled.not.registered.genes,
              'D.scaled.onlyR'=D.scaled.registered.genes,
              'D.registered.onlyR'=D.registered.registered.genes))
}

# sample.id.cols <- c('accession','delta_time')
# gene.col <- c('locus_name')
# expression.col <- 'mean_cpm'
# dt <- mean_df
#' @export
reformat_for_distance_calculation <- function(dt, sample.id.cols, gene.col, expression.col) {

  # concatenate sample.id columns to generate sample ids
  dt$sample.id <- dt[[sample.id.cols[1]]]
  if (length(sample.id.cols) > 1) {
    for (i in 2:length(sample.id.cols)) {
      # pad timepoint to 2 figures to help ordering
      if (class(dt[[sample.id.cols[i]]]) %in%  c('integer', 'numeric')) {
        dt[[sample.id.cols[i]]] <- stringr::str_pad(dt[[sample.id.cols[i]]], 2, pad='0')
      }

      dt$sample.id <- paste0(dt[['sample.id']], '-', dt[[sample.id.cols[i]]])
    }
  }

  # subset to just the relevant columns
  dt <- subset(dt, select=c('sample.id', gene.col, expression.col))

  # convert to wide format
  dt.w <- data.table::dcast(dt, locus_name ~ sample.id, value.var=eval(expression.col))

  return(dt.w)
}

#' @export
calc_sample_distance <- function(dt, gene.col) {
  # wrapper for calculate_pairwise_sample_distance
  # to calculate distance for all compared samples

  data.cols <- names(dt)[names(dt)!=eval(gene.col)]
  print(data.cols)

  i.cols <- c()
  j.cols <- c()
  ds <- c()

  for (i in 1:(length(data.cols)-1)) {
    i.col <- data.cols[i]
    for (j in (i+1):length(data.cols)) {
      j.col <- data.cols[j]

      curr.dt <- subset(dt, select=c(i.col, j.col))
      d <- calculate_pairwise_sample_distance_simple(curr.dt)

      i.cols <- c(i.cols, i.col, j.col) # distance metric is symmetrical
      j.cols <- c(j.cols, j.col, i.col)
      ds <- c(ds, d, d)
    }
  }
  out.df <- data.table::data.table(data.frame('x.sample'=i.cols, 'y.sample'=j.cols, 'distance'=ds))
}


#' Calculate pairwise sample distance
#'
#' @param df Dataframe contains the wide format expression of the two samples to be compared. Only genes which have data in both samples are considered (only relevant for registration case).
#'
#' @return Dataframe contains the distance between two samples
#' @export
calculate_pairwise_sample_distance_main <- function(df) {

  # Omit NA in the data
  df <- stats::na.omit(df)

  # Calculate distance
  df$sq_diff <- (df[, 1] - df[, 2])^2

  d <- mean(df$sq_diff)

  return(d)

}
