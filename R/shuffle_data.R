# unrandom.mean_df <- mean_df
# unrandom.all.df <- all.data.df
#' @export
shuffle_ro18_timepoints <- function(unrandom.mean_df, unrandom.all.df) {
  # shuffle the timepoints for each ro18 gene

  # split the mean_df
  col.df <- unrandom.mean_df[unrandom.mean_df$accession=='Col0',]
  ro18.df <- unrandom.mean_df[unrandom.mean_df$accession=='Ro18',]
  shuffled.ro18.df <-copy(ro18.df)
  # split the all.df
  col.all.df <- unrandom.all.df[unrandom.all.df$accession=='Col0',]
  ro18.all.df <- unrandom.all.df[unrandom.all.df$accession=='Ro18',]
  shuffled.ro18.all.df <-copy(ro18.all.df)

  # for each gene, replace the timepoints with the same shuffled timepoints
  # in both dfs
  ro18.times <- unique(ro18.df$timepoint)
  curr.locus <- unique(ro18.df$locus_name)[1]
  for (curr.locus in unique(ro18.df$locus_name)) {
    # generate common shuffle times lookup for this curr.locus
    shuffle.times <- sample(ro18.times)

    mean.replacement.times <- sapply(ro18.df$timepoint[ro18.df$locus_name==curr.locus],
                                     function(x) shuffle.times[match(x, ro18.times)])
    shuffled.ro18.df$timepoint[shuffled.ro18.df$locus_name==curr.locus] <- mean.replacement.times

    all.replacement.times <- sapply(ro18.all.df$timepoint[ro18.all.df$locus_name==curr.locus],
                                    function(x) shuffle.times[match(x, ro18.times)])
    shuffled.ro18.all.df$timepoint[shuffled.ro18.all.df$locus_name==curr.locus] <- all.replacement.times
  }

  mean_df <- rbind(col.df, shuffled.ro18.df)
  all.df <- rbind(col.all.df, shuffled.ro18.all.df)
  return(list(mean_df, all.df))
}

#' @export
shuffle_ro18_gene_names <- function(mean_df, out.all.df) {
  # shuffle the identities of the genes in the brassica
  # can't just do shuffle, becuase need to preserve which timepoints are from the same gene
  out.mean_df <- data.table::copy(mean_df)
  out.all.df <- data.table::copy(out.all.df)

  # make the gene lookup table for the same shuffled genes for both
  brassica.genes <- unique(out.mean_df$locus_name[out.mean_df$accession=='Ro18'])
  shuffled.genes <- sample(brassica.genes)
  shuffle.gene.lookup <- data.table::data.table(data.frame('gene.id'=brassica.genes, 'shuffled.id'=shuffled.genes))

  # change the gene names for the mean_df
  out.mean_df <- swap_gene_names(out.mean_df, shuffle.gene.lookup)
  # change the gene names for the all.df
  out.all.df <- swap_gene_names(out.all.df, shuffle.gene.lookup)

  return(list(out.mean_df, out.all.df))
}

#' @export
swap_gene_names <- function(df, shuffle.gene.lookup) {
  replacement.genes <- sapply(df$locus_name[df$accession=='Ro18'],
                              function(x) shuffle.gene.lookup$shuffled.id[match(x, shuffle.gene.lookup$gene.id)])

  replacement.genes <- as.character(replacement.genes) # otherwise returns a factor or strings,
  # depending on versions
  df$locus_name[df$accession=='Ro18'] <- replacement.genes

  return(df)
}
