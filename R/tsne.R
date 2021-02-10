# DATA <- exp
# rds_file
# plot_name <- 'test'
#do_tSNE(shift_5, rds_file, paste0('filt_model_shifted_&_scaled', curr_shift), sumBrassicas)
#DATA <- shift_5
#plot_name <- paste0('filt_model_shifted_&_scaled', curr_shift)
# DATA <- exp
# rds_file
# plot_name <- 'model_filtered_genes_no_scale'
# sumBrassicas
do_tSNE <- function(DATA, rds_file, plot_name, sumBrassicas) {

  if (sumBrassicas==T) {
    save.dir <- paste0('../graphs/', rds_file, '/tSNE_plots_summed_bra_copies')
  } else {
    save.dir <- paste0('../graphs/', rds_file, '/tSNE_plots_seperate_bra_copies')
  }

  D <- DATA[, c('group', 'locus_name', 'norm.cpm')]
  names(D) <- c('group', 'locus_name', 'expression')

  D <- data.table::dcast(D, group~locus_name, value.var='expression')

  # get rid of any genes (columns) with any NA
  D <- data.frame(D)
  #D <- D[, colSums(is.na(D)) == 0]
  tmp <- D[, colSums(is.na(D)) != 0]
  if (ncol(tmp)!=0) {
    print('got NAs in D. Something to do with using symbols as ids - if multi symbols for same At gene, somehow can end up with different ones in the different species..?')
    stopifnot(1==2)
  }

  #tmp <- cbind(Labels, tmp)
  ### train the model
  Labels <- D[, 1]
  #M <- stats::na.omit(as.matrix(D[, -1]))
  M <- as.matrix(D[, -1])
  nrow(M)
  ncol(M)
  ###rescale the data in each column - otherwise bigger expressed genes dominate
  sc.M <- scale(M)
  ###TMP FOR TESTING
  #sc.M <- sc.M[, c(1:800)]
  ###END TMP FOR TESTING

  # check that we get mean of 0 and sd of 1
  colMeans(sc.M)  # faster version of apply(scaled.dat, 2, mean)
  apply(sc.M, 2, stats::sd)

  #perp <- 1

  # ensure graph file exists
  if (!dir.exists(save.dir)) {
    dir.create(save.dir)
  }

  perp <- 1
  for (perp in c(1, 2, 3, 4, 5, 6, 8)) {
    # set to 8000 iterations -≥ pretty consistently minimises the error to ~0.35
    tsne <- Rtsne::Rtsne(sc.M, dims=2, perplexity=perp, verbose=T, max_iter=16000, pca=F)
    p1 <- tSNEPlot(tsne, Labels) + ggplot2::guides(color=F)
    tsne <- Rtsne::Rtsne(sc.M, dims=2, perplexity=perp, verbose=T, max_iter=16000, pca=F)
    p2 <- tSNEPlot(tsne, Labels) + ggplot2::guides(color=F)
    tsne <- Rtsne::Rtsne(sc.M, dims=2, perplexity=perp, verbose=T, max_iter=16000, pca=F)
    p3 <- tSNEPlot(tsne, Labels) + ggplot2::guides(color=F)
    tsne <- Rtsne::Rtsne(sc.M, dims=2, perplexity=perp, verbose=T, max_iter=16000, pca=F)
    p4 <- tSNEPlot(tsne, Labels) + ggplot2::guides(color=F)
    p <- ggpubr::ggarrange(p1, p2, p3, p4,
                           labels=c('a', 'b', 'c', 'd'),
                           ncol=2, nrow=2)
    p
    ggplot2::ggsave(paste0(save.dir, '/', plot_name, '_p=', perp,'.pdf'), scale=1.5)
  }

  # try with PCA before tSNE, as it seems to give an ok gradient timecourse.
  # for (perp in c(2,3,4,5,68,12,15)) {
  #   # set to 8000 iterations -≥ pretty consistently minimises the error to ~0.35
  #   tsne <- Rtsne::Rtsne(sc.M, dims=2, perplexity=perp, verbose=T, max_iter=16000, pca=T, initial_dims=12)
  #   p1 <- tSNEPlot(tsne, Labels) + ggplot2::guides(color=F)
  #   tsne <- Rtsne::Rtsne(sc.M, dims=2, perplexity=perp, verbose=T, max_iter=16000, pca=T, initial_dims=12)
  #   p2 <- tSNEPlot(tsne, Labels) + ggplot2::guides(color=F)
  #   tsne <- Rtsne::Rtsne(sc.M, dims=2, perplexity=perp, verbose=T, max_iter=16000, pca=T, initial_dims=12)
  #   p3 <- tSNEPlot(tsne, Labels) + ggplot2::guides(color=F)
  #   tsne <- Rtsne::Rtsne(sc.M, dims=2, perplexity=perp, verbose=T, max_iter=16000, pca=T, initial_dims=12)
  #   p4 <- tSNEPlot(tsne, Labels) + ggplot2::guides(color=F)
  #   p <- ggpubr::ggarrange(p1, p2, p3, p4,
  #                  labels=c('a', 'b', 'c', 'd'),
  #                  ncol=2, nrow=2)
  #   ggplot2::ggsave(paste0('../graphs/', rds_file, '/tSNE_plots/', plot_name, '_p=', perp,'_PCA.pdf'), scale=1.5)
  # }
}

tSNEPlot <- function(tsne, Labels) {
  position <- data.frame(tsne$Y)
  position <- cbind(position, Labels)
  names(position) <- c('x', 'y', 'label')

  position <- data.table::data.table(position)
  position[, c('accession', 'timepoint', 'letter'):=data.table::tstrsplit(label, split='-')]
  p <- ggplot2::ggplot(position)+
    ggplot2::aes(x=x, y=y, color=accession, label=paste(timepoint, letter, sep='-')) +
    ggplot2::geom_point()+
    ggrepel::geom_text_repel(size=2)
  return(p)
}
