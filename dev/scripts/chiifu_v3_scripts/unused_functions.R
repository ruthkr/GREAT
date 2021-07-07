# Unused functions ----

# real.model.comparisons <- real.registered
# shuffled.model.comparisons <- shuffled.registered
make_barplot_of_shifts <- function(real.model.comparisons, shuffled.model.comparisons) {
  # plot the shifts applied in real and shuffled data.
  # NB THIS DOES NOT CUT DOWN TO ONLY INCLUDE ONES WHERE REGISTRATION IS THE BETTER MODEL!

  # do the calculations for the shuffled.shifts uncertainty in counts
  shuffled.model.comparisons[, count:=.N, by=.(stretch, shift, job)]

  shuffled.reg.counts <- unique(shuffled.model.comparisons[, c('stretch', 'shift','job', 'count')])


  shuffled.4.plot <- shuffled.reg.counts[, .('mean.count'=mean(count),'sd.count'=sd(count)), by=.(stretch, shift)]
  shuffled.4.plot$is.real <- 'shuffled homologue mapping'

  real.4.plot <- real.model.comparisons[, .('mean.count'=.N, 'sd.count'=NA), by=.(stretch, shift)]
  real.4.plot$is.real <- 'real homologue mapping'

  both.4.plot <- rbind(real.4.plot, shuffled.4.plot)
  both.4.plot$stretch <- paste0('stretch = ', both.4.plot$stretch, 'x')
  both.4.plot$stretch <- factor(both.4.plot$stretch, levels=c('stretch = 1x', 'stretch = 1.5x', 'stretch = 2x'))

  p <- ggplot(both.4.plot, aes(x=shift, y=mean.count, fill=is.real, color=is.real))+
    geom_bar(stat='identity', alpha=0.5, position=position_dodge(width=0.4))+
    geom_errorbar(aes(ymin=mean.count-2*sd.count, ymax=mean.count+2*sd.count),
                  width=0.2, position=position_dodge(width=0.4))+
    facet_wrap(~stretch, ncol=1)+
    ylab('count')+
    scale_x_continuous(breaks=scales::pretty_breaks())+
    theme_bw()+
    guides(fill=guide_legend(nrow=2,byrow=TRUE),
           color=guide_legend(nrow=2,byrow=TRUE))+
    theme(legend.position='top',
          axis.title = element_text(size=10),
          axis.text = element_text(size=6),
          strip.text = element_text(size=8),
          legend.title = element_blank(),
          legend.margin = margin(0,0,-9,-4))

  return(p)
}

#real.registered <- real.model.comparisons
make_barplot_of_shifts_real_only <- function(real.model.comparisons) {
  # plot the shifts applied in real and shuffled data.
  # NB THIS DOES NOT CUT DOWN TO ONLY INCLUDE ONES WHERE REGISTRATION IS THE BETTER MODEL!

  real.4.plot <- real.model.comparisons[, .('mean.count'=.N, 'sd.count'=NA), by=.(stretch, shift)]
  real.4.plot$is.real <- 'real homologue mapping'

  both.4.plot <- real.4.plot
  both.4.plot$stretch <- paste0('stretch = ', both.4.plot$stretch, 'x')
  both.4.plot$stretch <- factor(both.4.plot$stretch, levels=c('stretch = 1x', 'stretch = 1.5x', 'stretch = 2x'))

  p <- ggplot(both.4.plot, aes(x=shift, y=mean.count))+
    geom_segment(aes(x=shift, xend=shift, y=0, yend=mean.count), size=4)+
    facet_wrap(~stretch, ncol=1)+
    ylab('count')+
    xlab('shift (d)')+
    #scale_x_continuous(breaks=scales::pretty_breaks())+
    theme_bw()+
    theme(legend.position='none',
          axis.title = element_text(size=10),
          axis.text = element_text(size=6),
          strip.text = element_text(size=8),
          legend.title = element_blank(),
          legend.margin = margin(0,0,-9,-4))

  return(p)
}

calculate_between_sample_distance_for_shuffled_data <- function(shuffled.data.dir) {
  # wrapper for calculate_between_sample_distance() which applies it to all the shuffled results

  # print('loading shuffled data...')
  # all.mean_df <- load_shuffled_data(shuffled.data.dir, 'mean_df')
  # all.sc.df <- load_shuffled_data(shuffled.data.dir, 'mean.sc')
  # all.registered.df <- load_shuffled_data(shuffled.data.dir, 'imputed.mean_df')

  print('calculating distances...')
  jobIds <- get_job_suffixes(shuffled.data.dir)

  D.mean.list <- rep(list(0), length(jobIds))
  D.sc.list <- rep(list(0), length(jobIds))
  D.registered.list <- rep(list(0), length(jobIds))
  D.scaled.onlyNR.list <- rep(list(0), length(jobIds))
  D.scaled.onlyR.list <- rep(list(0), length(jobIds))
  D.registered.onlyR.list <- rep(list(0), length(jobIds))
  for (i in 1:length(jobIds)) {
    if (i %% 5==0 ){
      print(paste0(i, ' / ', length(jobIds)))
    }
    curr.job <- jobIds[i]

    curr.mean_df <- readRDS(paste0(shuffled.data.dir, 'mean_df_', curr.job, '.rds')) #all.mean_df[all.mean_df$job==curr.job,]
    curr.sc.df <- readRDS(paste0(shuffled.data.dir, 'mean_df.sc_', curr.job, '.rds')) #all.sc.df[all.sc.df$job==curr.job,]
    curr.imputed.mean_df <- readRDS(paste0(shuffled.data.dir, 'imputed.mean_df_', curr.job, '.rds')) #all.registered.df[all.registered.df$job==curr.job,]

    O <- calculate_between_sample_distance(curr.mean_df, curr.sc.df, curr.imputed.mean_df)
    curr.D.mean <- O[['D.mean']]
    curr.D.scaled <- O[['D.scaled']]
    curr.D.registered <- O[['D.registered']]
    curr.D.scaled.onlyNR <- O[['D.scaled.onlyNR']]
    curr.D.scaled.onlyR <- O[['D.scaled.onlyR']]
    curr.D.registered.onlyR <- O[['D.registered.onlyR']]

    curr.D.mean$title <- 'mean expression'
    curr.D.scaled$title <- 'scaled expression (all)'
    curr.D.registered$title <- 'registered expression (all)'
    curr.D.scaled.onlyNR$title <- 'scaled expression (not registered only)'
    curr.D.scaled.onlyR$title <- 'scaled expression (registered only)'
    curr.D.registered.onlyR$title <- 'registered expression (registered only)'

    curr.D.mean$ID <- curr.job
    curr.D.scaled$ID <- curr.job
    curr.D.registered$ID <- curr.job
    curr.D.scaled.onlyNR$ID <- curr.job
    curr.D.scaled.onlyR$ID <- curr.job
    curr.D.registered.onlyR$ID <- curr.job

    D.mean.list[[i]] <- curr.D.mean
    D.sc.list[[i]] <- curr.D.scaled
    D.registered.list[[i]] <- curr.D.registered
    D.scaled.onlyNR.list[[i]] <- curr.D.scaled.onlyNR
    D.scaled.onlyR.list[[i]] <- curr.D.scaled.onlyR
    D.registered.onlyR.list[[i]] <- curr.D.registered.onlyR
  }
  # package the list of Ds to dfs for each condition
  print('rbinding lists to dt...')
  D.mean.shuffled <- do.call('rbind', D.mean.list)
  D.scaled.shuffled <- do.call('rbind', D.sc.list)
  D.registered.shuffled <- do.call('rbind', D.registered.list)
  D.scaled.onlyNR.shuffled <- do.call('rbind', D.scaled.onlyNR.list)
  D.scaled.onlyR.shuffled <- do.call('rbind', D.scaled.onlyR.list)
  D.registered.onlyR.shuffled <- do.call('rbind', D.registered.onlyR.list)

  O <- list('D.mean.shuffled'=D.mean.shuffled,
            'D.scaled.shuffled'=D.scaled.shuffled,
            'D.registered.shuffled'=D.registered.shuffled,
            'D.scaled.onlyNR.shuffled'=D.scaled.onlyNR.shuffled,
            'D.scaled.onlyR.shuffled'=D.scaled.onlyR.shuffled,
            'D.registered.onlyR.shuffled'=D.registered.onlyR.shuffled
  )
  return(O)
}

plot_number_of_registered_genes_in_shuffled <- function(real.data.dir, shuffled.data.dir) {
  # plot the number of registered genes in the real data, and in the shuffled
  # data. Plot the number in the real data as a vertical line, and calculate the
  # quantile.

  # get real number of genes registered
  num.registered.real <- get_num_registered_genes(paste0(real.data.dir, 'imputed.mean_df.rds'))
  # get vector of number of genes registered after random shuffling
  imputed.random.files <- list.files(path=shuffled.data.dir, pattern='imputed')
  imputed.random.paths <- paste0(shuffled.data.dir, imputed.random.files)
  num.registered.random <- sapply(imputed.random.paths, get_num_registered_genes)
  num.registered.random <- data.frame(as.numeric(num.registered.random))
  names(num.registered.random) <- 'num.registered'
  # calculate quantile
  real.quantile <- ecdf_fun(num.registered.random$num.registered, num.registered.real)
  # make plot

  Z.score <- (num.registered.real - mean(num.registered.random$num.registered)) / sd(num.registered.random$num.registered)
  pval <- pnorm(-abs(Z.score))
  print('P value, assuming normal distribution:')
  print(pval)

  p <- ggplot(num.registered.random, aes(x=num.registered))+
    geom_histogram(aes(y=..density..), bins=30, alpha=0.8)+
    geom_vline(aes(xintercept=num.registered.real), size=1.2)+
    stat_function(fun=dnorm, args=list(mean=mean(num.registered.random$num.registered),
                                       sd=sd(num.registered.random$num.registered)),
                  size=0.5)+
    theme_bw()+
    xlab(paste0('number of genes identified as having similar gene\nexpression profiles in both organisms'))+
    theme(
      axis.title = element_text(size=10),
      axis.text = element_text(size=6)
    )
  #ggtitle(paste0('Quantile : ', real.quantile))
  return(p)
}

plot_goI_expression <- function(summed.GoIs.df) {

  # make plot of gene expression in Col0, and in Ro18.
  # truncate data so can see expression nicely on the same scale
  summed.GoIs.df <- summed.GoIs.df[summed.GoIs.df$timepoint <=21,]

  summed.GoIs.df[, scaled.cpm:=my_scale(mean_cpm), by=.(Ara.id, accession)]
  morphology.equiv.df <- data.frame('accession'=c('Col-0', 'R-O-18'), 'floral.transition.time'=c(14, 35))
  summed.GoIs.df$accession <- as.character(summed.GoIs.df$accession)
  summed.GoIs.df$accession[summed.GoIs.df$accession=='Col0'] <- 'Col-0'
  summed.GoIs.df$accession[summed.GoIs.df$accession=='Ro18'] <- 'R-O-18'

  curr.acc <- 'Col-0'
  plot.list <- list()
  for (curr.acc in c('Col-0', 'R-O-18')) {

    curr.p <- ggplot(summed.GoIs.df[summed.GoIs.df$accession==curr.acc,],
                     aes(x=timepoint, y=mean_cpm, color=Ara.name, fill=Ara.name))+
      #geom_vline(data=summed.GoIs.df[summed.GoIs.df$accession=='Col-0' & summed.GoIs.df$accession==curr.acc,], aes(xintercept=floral.transition.time), size=1)+
      #geom_vline(data=summed.GoIs.df[summed.GoIs.df$accession=='R-O-18' & summed.GoIs.df$accession==curr.acc,,], aes(xintercept=floral.transition.time), size=1)+
      geom_vline(data=morphology.equiv.df[morphology.equiv.df==curr.acc,], aes(xintercept=floral.transition.time), size=1)+

      stat_summary(fun=mean, geom='line', size=1)+
      stat_summary(fun.data=mean_se, fun.args=list(mult=1.96),geom='ribbon',
                   color=NA, alpha=0.3)+
      geom_point(size=0.4)+
      facet_wrap(~accession, ncol=1, scales='free')+
      scale_x_continuous(breaks=scales::pretty_breaks())+
      theme_bw()+
      ylab('expression (cpm)')+
      xlab('timepoint (d)')+
      guides(fill=guide_legend(nrow=2,byrow=TRUE),
             color=guide_legend(nrow=2,byrow=TRUE))+
      theme(legend.position='top',
            legend.title = element_blank(),
            legend.text = element_text(face='italic', size=8),
            axis.title = element_text(size=10),
            axis.text=element_text(size=6),
            legend.margin=margin(0,0,-10,0))

    leg <- get_legend(curr.p)
    leg.plot <- as_ggplot(leg)
    leg.plot + theme(plot.margin=c(0,0,0,0))

    curr.p <- curr.p + theme(legend.position='none')

    if (curr.acc == 'Col-0') {
      curr.p <- curr.p + theme(axis.title.x=element_blank())+
        coord_cartesian(ylim=c(0, 250))
    }

    plot.list <- c(plot.list, list(curr.p))
  }
  plot.list <- c(list(leg.plot), plot.list)
  p.grid <- plot_grid(plotlist=plot.list, rel_heights = c(0.3,0.93,1), ncol=1, align='v')

  return(p.grid)
}

plot_registered_GoIs_for_comparible_timepoints <- function(all.stretched.df) {

  # make plot of gene expression after registration - only plot compared timepoints
  # cut down to compared timepoints for each gene

  registered.plot.df <- all.stretched.df
  AGL24.df <- registered.plot.df[registered.plot.df$locus_name=='AGL24' &
                                   (registered.plot.df$shifted_time <=35) &
                                   (registered.plot.df$shifted_time >=13),]
  AP1.df <- registered.plot.df[registered.plot.df$locus_name=='AP1' &
                                 (registered.plot.df$shifted_time <=25) &
                                 (registered.plot.df$shifted_time >=11),]
  AP3.df <- registered.plot.df[registered.plot.df$locus_name=='AP3' &
                                 (registered.plot.df$shifted_time <=25) &
                                 (registered.plot.df$shifted_time >=13),]
  LFY.df <- registered.plot.df[registered.plot.df$locus_name=='LFY' &
                                 (registered.plot.df$shifted_time <=19) &
                                 (registered.plot.df$shifted_time >=10),]
  SOC1.df <- registered.plot.df[registered.plot.df$locus_name=='SOC1' &
                                  (registered.plot.df$shifted_time <=31) &
                                  (registered.plot.df$shifted_time >=10),]
  registered.plot.df <- rbind(AGL24.df, AP1.df, AP3.df, LFY.df, SOC1.df)

  registered.plot.df$accession <- as.character(registered.plot.df$accession)
  registered.plot.df$accession[registered.plot.df$accession=='Col0'] <- 'Col-0'
  registered.plot.df$accession[registered.plot.df$accession=='Ro18'] <- 'DH'

  p.registered <- ggplot(registered.plot.df, aes(x=shifted_time, y=mean_cpm, color=accession, fill=accession))+
    stat_summary(fun=mean, geom='line', size=1)+
    stat_summary(fun.data=mean_se, fun.args=list(mult=1.96),geom='ribbon',
                 color=NA, alpha=0.3)+
    geom_point(size=0.4)+
    theme_bw()+
    xlab('registered time (d)')+
    ylab('normalised expression')+
    facet_wrap(~locus_name, scales='free', ncol=2)+
    scale_x_continuous(breaks=scales::pretty_breaks())+
    theme(legend.position = 'top',
          legend.title = element_blank(),
          axis.title = element_text(size=10),
          axis.text=element_text(size=6),
          strip.text=element_text(face='italic'),
          legend.margin=margin(21,0,0,0))

  return(p.registered)
}

make_quantile_heatmaps <- function(Q.mean, Q.scaled, Q.registered) {

  p.mean <- make_heatmap_quantile(Q.mean, 'mean expression')
  p.scaled <- make_heatmap_quantile(Q.scaled, 'scaled mean expression')
  p.registered <- make_heatmap_quantile(Q.registered, 'registered & scaled mean expression')

  p.all <- plot_grid(p.mean, p.scaled, p.registered, ncol=1)

  return(p.all)
}

ecdf_fun <- function(x,perc) {
  # estimate the quantile of the perc value in the x distribution
  # https://stats.stackexchange.com/questions/50080/estimate-quantile-of-value-in-a-vector/114493
  stats::ecdf(x)(perc)
}

# real.dt <- D.mean.real
# shuffled.dt <- D.mean.shuffled
get_all_quantiles <- function(real.dt, shuffled.dt) {
  # for each real pairwise distance, calculate the quantile of that
  # in the distribution of distances for that comparison in the randomly
  # shuffled gene expression runs.

  out.dt <- data.table(data.frame('x.sample'=real.dt$x.sample,
                                  'y.sample'=real.dt$y.sample,
                                  'quantile'=NA,
                                  'title'=real.dt$title))

  # x.val <- unique(real.dt$x.sample)[1]
  # y.val <- unique(real.dt$y.sample)[1]
  for (x.val in unique(real.dt$x.sample)) {
    for (y.val in unique(real.dt$y.sample)) {
      if ( (x.val == y.val) |
           (grepl('Col0', x.val) & grepl('Col0', y.val)) |
           (grepl('Ro18', x.val) & grepl('Ro18', y.val)) ) {
        quantile <- NA
      } else {
        # real distance for this comparison
        real.dist <- real.dt$distance[real.dt$x.sample==x.val &
                                        real.dt$y.sample==y.val]

        # distribution of shuffled distances for this comparison
        shuffled.distances <- shuffled.dt$distance[shuffled.dt$x.sample==x.val &
                                                     shuffled.dt$y.sample==y.val]

        if (is.na(real.dist) | length(unique(shuffled.distances))==1) {
          quantile <- NA
        } else {
          #print(real.dist)
          #print(length(shuffled.distances))
          quantile <- ecdf_fun(shuffled.distances, real.dist)
        }
      }

      out.dt$quantile[out.dt$x.sample==x.val &
                        out.dt$y.sample==y.val] <- quantile
    }
  }
  return(out.dt)
}

make_shuffled_data_heatmaps <- function(D.mean, D.scaled, D.registered) {

  p.mean <- make_heatmap_w_shuffled(D.mean, 'mean expression')
  p.scaled <- make_heatmap_w_shuffled(D.scaled, 'scaled mean expression')
  p.registered <- make_heatmap_w_shuffled(D.registered, 'registered & scaled mean expression')

  p.all <- plot_grid(p.mean, p.scaled, p.registered, ncol=1)

  return(p.all)
}

load_shuffled_distances <- function(dir.name, D.type, numJobs, numRuns) {
  # dir.name the directory to look in (e.g. "rescaled_shuffled_distance") in intermediate_data/gene_registration
  # D.type in c('D.mean', 'D.scaled', 'D.registered') depending
  # how the shuffled expression data was transformed
  # mean of biological replicates
  # that + expression scaled(center=T, scale=T)
  # that + registered (scale using compared points only)

  # numJobs : how many separate jobs were submitted
  # numRuns : how many runs each job did

  file.stem <- paste0('../intermediate_data/gene_registration/', dir.name, '/', D.type)

  listOfDfs <- vector("list", length=numJobs * numRuns)
  idx <- 1
  for (jobID in 1:numJobs) {
    for (runID in 1:numRuns) {
      curr.path <- paste0(file.stem, '_', jobID, '_', runID, '.rds')
      #print(curr.path)
      curr.df <- readRDS(curr.path)
      curr.df$ID <- paste0(jobID, '_', runID)

      listOfDfs[[idx]] <- curr.df
      idx <- idx + 1
    }
  }
  Df <- do.call('rbind', listOfDfs)

  return(Df)
}


# ro18_rds_file <- '../final_data/rds/ro18_leaf_reannotated.rds'
# sari14_rds_file <- '../final_data/rds/sari14_chiifu_leaf_Jul2020.rds'
load_mean_df_sari_ro18 <- function(ro18_rds_file, sari_rds_file, pCor.th) {


  ro18.exp <- readRDS(ro18_rds_file)
  #ro18.exp$FPKM <- NULL
  sari14.exp <- readRDS(sari14_rds_file)
  #sari14.exp$TPM <- NULL

  # check that comparing the same tissue. Can't easily check
  # that the same alignment though.
  if (unique(ro18.exp$tissue != unique(sari14.exp$tissue))) {
    print("compared tissues aren't the same in sari and exp")
    stop()
  }

  # get the mean of each gene at each timepoint
  ro18.exp[, mean_cpm:=mean(norm.cpm), by=.(CDS.model, accession, tissue, timepoint)]
  sari14.exp[, mean_cpm:=mean(norm.cpm), by=.(CDS.model, accession, tissue, timepoint)]

  exp <- rbind(ro18.exp, sari14.exp)
  mean_df <- unique(exp[, c('CDS.model', 'accession', 'tissue', 'timepoint', 'mean_cpm')])
  names(mean_df) <-  c('locus_name', 'accession', 'tissue', 'timepoint', 'mean_cpm')
  names(exp)[names(exp)=='CDS.model'] <- 'locus_name'
  # checking not cutting out key genes - need to know what they're called now
  # FT
  # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.8543'),]
  # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.41993'),]
  # # SOC1
  # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.15712'),]
  # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.26487'),]

  # ggplot(tmp, aes(x=timepoint, y=mean_cpm))+
  #   geom_line()

  # filter mean_df to remove genes with very low expression - remove if max is less than 5, and less than half timepoints expressed
  # greater than 1 in both accessions
  ro18.df <- mean_df[mean_df$accession=='Ro18']
  ro18.df[, keep:=(max(mean_cpm) > 2 | mean(mean_cpm > 1) > 0.5) , by=.(locus_name)]
  ro18.keep.genes <- unique(ro18.df$locus_name[ro18.df$keep==TRUE])
  #'MSTRG.8543' %in% ro18.keep.genes

  sari14.df <- mean_df[mean_df$accession=='sarisha14']
  sari14.df[, keep:=(max(mean_cpm) > 2 | mean(mean_cpm > 1) > 0.5) , by=.(locus_name)]
  sari14.keep.genes <- unique(sari14.df$locus_name[sari14.df$keep==TRUE])
  #'MSTRG.8543' %in% sari14.keep.genes

  keep.genes <- intersect(ro18.keep.genes, sari14.keep.genes)
  mean_df <- mean_df[mean_df$locus_name %in% keep.genes,]


  # filter to remove genes with low correlation between individual timepoints
  # and mean timepoints in either accession
  keep.genes <- filter.low.variability(exp, pCor.th)
  #'MSTRG.8543' %in% keep.genes

  mean_df <- mean_df[mean_df$locus_name %in% keep.genes, ]
  #'MSTRG.8543' %in% mean_df$locus_name

  # # checking not cutting out key genes - need to know what they're called now
  # # FT
  # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.8543'),]
  # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.41993'),]
  # # SOC1
  # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.15712'),]
  # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.26487'),]
  #

  print(paste0(length(unique(mean_df$locus_name)), ' genes considered in the comparison'))
  rm(ro18.df, keep.genes)

  exp <- exp[exp$locus_name %in% unique(mean_df$locus_name)]
  exp <- subset(exp, select=c('locus_name', 'accession', 'tissue', 'timepoint',
                              'norm.cpm', 'group'))
  names(exp)[names(exp)=='norm.cpm'] <- 'mean_cpm'
  return(list(mean_df, exp))
}

plot_all_real_distance_heatmaps <- function(D.mean, D.scaled, D.scaled.onlyNR, D.scaled.onlyR, D.registered, D.registered.onlyR) {

  # change title used for facetting
  D.mean$title <- 'measured expression'
  D.scaled$title <- 'scaled expression'
  D.scaled.onlyNR$title <- 'only not registered'
  D.scaled.onlyR$title <- 'only registered'
  D.registered$title <- 'registered expression'
  D.registered.onlyR$title <- 'only registered'


  # MAKE THE BITS of the PLOTS FOR THE REAL DATA
  p.blank <- ggplot() + theme_void() # blank facet used for layout

  # for RAW DATA, show all quadrants
  p.mean <- make_heatmap(D.mean, '')
  #p.mean.grid <- plot_grid(p.blank, p.mean, p.blank, nrow=1, rel_widths = c(0.1, 2,0.1))
  p.mean.grid <- plot_grid(p.mean, labels=c('a'), label_size = 18, label_fontface = 'plain')

  # for SCALED DATA, show only lower right
  # cut down data to only show lower right quadrant
  D.scaled.all <- D.scaled[grepl(pattern='Col0', D.scaled$y.sample) &
                             grepl(pattern='Ro18', D.scaled$x.sample),]
  D.scaled.onlyNR <- D.scaled.onlyNR[grepl(pattern='Col0', D.scaled.onlyNR$y.sample) &
                                       grepl(pattern='Ro18', D.scaled.onlyNR$x.sample),]
  D.scaled.onlyR <- D.scaled.onlyR[grepl(pattern='Col0', D.scaled.onlyR$y.sample) &
                                     grepl(pattern='Ro18', D.scaled.onlyR$x.sample),]
  D.scaled <- rbind(D.scaled.all)#, D.scaled.onlyNR, D.scaled.onlyR)
  #D.scaled$title <- factor(D.scaled$title, 'all genes')#, 'only registered', 'only not registered'))

  p.scaled.grid <- make_heatmap(D.scaled.all, '')

  # for REGISTERED DATA
  D.registered.all <- D.registered[grepl(pattern='Col0', D.registered$y.sample) &
                                     grepl(pattern='Ro18', D.registered$x.sample),]
  D.registered.onlyR <- D.registered.onlyR[grepl(pattern='Col0', D.registered.onlyR$y.sample) &
                                             grepl(pattern='Ro18', D.registered.onlyR$x.sample),]
  D.registered <- rbind(D.registered.all)#, D.registered.onlyR)
  D.registered <- D.registered[!(D.registered$y.sample %in% c('Col0-34', 'Col0-35'))]

  p.registered.g <- make_heatmap(D.registered, '', y.axis.fontsize=6)
  #p.registered.grid <- plot_grid(p.registered.g, p.blank, nrow=1, rel_widths = c(0.7,0.3))

  bcgrid <- plot_grid(p.scaled.grid, p.registered.g, nrow=1, labels=c('b', 'c'), label_size = 18, label_fontface = 'plain')
  bcgrid



  # combine with no treatment
  p.grid <- plot_grid(p.mean.grid, bcgrid, ncol=1, align='v', rel_heights = c(0.8, 0.5))
  p.grid
  return(p.grid)
}

plot_registration_for_exemplar_genes <- function(all.rep.shifted.data, GoIs) {
  # cut down to genes to plot
  tmp <- all.rep.shifted.data[all.rep.shifted.data$locus_name %in% GoIs]

  p.separate <- ggplot(tmp, aes(x=timepoint, y=mean_cpm, color=locus_name, fill=locus_name,
                                shape=accession, linetype=accession))+
    stat_summary(fun=mean, geom='line', size=1)+
    #stat_summary(fun=mean, geom='point')+
    stat_summary(fun.data=mean_se, fun.args=list(mult=1.96),geom='ribbon',
                 color=NA, alpha=0.3)+
    geom_point(size=0.5)+
    facet_wrap(~accession, ncol=1, scales='free_x')+
    theme_bw()+
    ylab('scaled gene expression')+
    xlab('timepoint (d)')+
    guides(fill=guide_legend(title='Gene', title.position='top',
                             override.aes=list(size=1)),
           color=guide_legend(title='Gene', title.position='top',
                              override.aes=list(size=1)),
           shape=FALSE, linetype=FALSE)+
    theme(legend.position='top',
          legend.title=element_text(size=10),
          legend.text=element_text(size=6),
          legend.margin=margin(0,0,0,0),
          legend.box.margin = margin(0,0,-10,-25),
          legend.spacing.x = unit(0.01, 'cm'),
          axis.text = element_text(size=6),
          axis.title = element_text(size=8),
          strip.text = element_text(size=10),
          strip.background = element_rect(fill='white'))


  tmp$label <- 'Registered'
  p.shifted <- ggplot(tmp, aes(x=shifted_time, y=mean_cpm,
                               color=locus_name, fill=locus_name,
                               shape=accession, linetype=accession))+
    stat_summary(fun=mean, geom='line', size=1)+
    #stat_summary(fun=mean, geom='point')+
    stat_summary(fun.data=mean_se, fun.args=list(mult=1.96),
                 geom='ribbon', color=NA, alpha=0.15)+
    geom_point(size=0.5)+
    facet_wrap(~label)+
    theme_bw()+
    ylab('scaled gene expression')+
    xlab('registered timepoint')+
    theme(legend.position='none',
          axis.text = element_text(size=6),
          axis.title = element_text(size=8),
          strip.text = element_text(size=10),
          strip.background = element_rect(fill='white'))

  p.grid <- plot_grid(p.separate, p.shifted, ncol=1, rel_heights = c(0.6, 0.4))
  # p.grid
  # ggsave('./test.pdf',
  #        plot=p.grid,
  #        width=3.2, height=6)

  return(p.grid)
}

#dt <- subset(mean.dt.sc.w, select=c('Col0-00', 'Col0-04'))
calculate_pairwise_sample_distance <- function(dt) {
  # dt is the wide format expression of the two samples to be compared.
  # only genes which have data in both samples are considered (only
  # relevant for registration case).

  # filter to genes with data in both.
  dt <- na.omit(dt)
  # calculate distance
  dt$sq.diff <- (dt[, 1] - dt[, 2])^2
  # calculate mean of each gene
  dt$mean <- (dt[, 1] + dt[, 2]) / 2
  dt$mean[dt$mean==0] <- 1 # divide by mean for score
  # calculate score for each gene
  dt$score <- dt$sq.diff / abs(dt$mean)

  d <- mean(dt$score)
  return(d)
}

rescale.by.mean <- function(v) {
  if (mean(v) == 0) {
    return(v)
  } else {
    return (v / mean(v))
  }
}

get_data_unique_symbol <- function(brassica_name) {

  rds_k <- 'klepikova' # 'klepikova

  rdsdir <- 'final_data/rds/'
  # rdsdir <- '../final_data/rds/'
  rdspath_k <- paste0(rdsdir, rds_k, '.rds')
  rdspath_b <- paste0(rdsdir, brassica_name, '.rds')

  expression_b <- readRDS(rdspath_b)
  expression_k <- readRDS(rdspath_k)

  # id data
  ID_TABLE <- readRDS(paste0('../reference_data/ID_TABLE_brapa-v3.rds'))
  ID_TABLE <- unique(ID_TABLE[, c('locus_name', 'symbol', 'CDS.model')])

  # this screws up the number of copies of each gene unless reorder based on symbol name first!
  # ID_TABLE <- ID_TABLE[order(ID_TABLE$symbol, decreasing=T), ]
  # ID_TABLE <- ID_TABLE[!duplicated(ID_TABLE[,c('locus_name', 'CDS.model')]), ]
  # length(unique(ID_TABLE$CDS.model))


  # add symbol info and cut down
  little.IDT <- ID_TABLE[, c('CDS.model', 'symbol')]
  expression_b <- merge(expression_b, little.IDT, by='CDS.model', allow.cartesian = T)

  little.IDT <- ID_TABLE[, c('locus_name', 'symbol')]
  names(little.IDT) <- c('CDS.model', 'symbol')
  expression_k <- unique(merge(expression_k, little.IDT, by='CDS.model', allow.cartesian = T))

  # cut down to only have genes with symbol present in both datasets
  common_symbols <- intersect(expression_k$symbol, expression_b$symbol)
  common_symbols <- common_symbols[!is.na(common_symbols)]
  common_symbols <- common_symbols[common_symbols!='']
  # get rid of symbols which are actually the same underlying gene model
  idt <- ID_TABLE[ID_TABLE$symbol %in% common_symbols]
  idt <- idt[order(idt$symbol, decreasing=T), ]
  idt <- idt[!duplicated(idt[, c('locus_name', 'CDS.model')]), ]
  common_symbols <- unique(idt$symbol)

  expression_k <- expression_k[expression_k$symbol %in% common_symbols, ]
  expression_b <- expression_b[expression_b$symbol %in% common_symbols, ]

  # join the two datasets into 1 & housekeeping
  expression <- rbind(expression_b, expression_k)

  # cut down to only have tissue present in both dataset
  #expression <- expression[expression$tissue=='apex', ]

  # cut down to only have similar/equivalent timepoint
  # cut down extra Col0 times
  #expression <- expression[expression$accession!='Col0' | expression$timepoint <=14, ]
  # cut down extra ZS11 times
  # expression <- expression[expression$accession!='ZS11' | (expression$timepoint<63 & expression$timepoint >=43), ]
  # # cut out dodgy 51 day timepoint
  # expression <- expression[expression$dataset!='ds_2018_06_19' | expression$timepoint != 51, ]

  rm(expression_b, expression_k, ID_TABLE, little.IDT)
  gc()

  return(expression)
}

same_dev_rescale <- function(x) {
  rs <- x$norm.cpm / mean(x$norm.cpm)
  # if (unique(x$accession)[1] != 'Col0') {
  #   rs <- x$norm.cpm / mean(x$norm.cpm)
  #   #rs <- scale(x)
  # } else if (unique(x$accession)[1] == 'Col0') {
  #   # rescale by the mean of the timepoints think are present in brassica too.
  #   k <- x$norm.cpm[x$timepoint <=17] # worked clearly with 13!
  #   rs <- x$norm.cpm / mean(k)
  # }
  return(rs)
}

calculate_sample_distance <- function(exp) {
  # want to calculate the average euclidean distance
  #exp$tmp.label <- paste0(exp$accession, '-', exp$timepoint)
  i <- unique(exp$group)[1]
  dist.list <- list()
  for (i in unique(exp$group)) {
    i.expression <- exp[group==i, c('locus_name', 'norm.cpm')]
    #i.mean.expression <- i.expression[, .(norm.cpm = mean(norm.cpm)), by=.(locus_name)]
    j <- unique(exp$group)[2]
    for (j in unique(exp$group)) {
      j.expression <- exp[group==j, c('locus_name', 'norm.cpm')]
      #j.mean.expression <- j.expression[, .(norm.cpm = mean(norm.cpm)), by=.(locus_name)]
      curr.expression <- merge(i.expression, j.expression, by='locus_name')
      curr.expression$diff <- sqrt((curr.expression$norm.cpm.x - curr.expression$norm.cpm.y)**2)
      dist <- mean(curr.expression$diff)
      out.df <- data.frame('x'=i, 'y'=j, 'distance'=dist)
      dist.list <- c(dist.list, list(out.df))
    }
  }
  dist.list[[1]]
  dist.df <- do.call('rbind', dist.list)
  return(dist.df)
}

#test <- mean_df
get_best_shift_old <- function(curr_sym, test) {
  # here statistic used to assess best is 1 / mean(squared difference)

  min_shift <- 5

  # cut to get a single symbol
  test <- test[test$locus_name==curr_sym, ]

  # get the mean expression vectors for the current gene
  AtVec <- test$mean_cpm[test$accession=='Col0']
  BrVec <- test$mean_cpm[test$accession=='Ro18']

  # get the arabidopsis, brassica timepoints will to compare to Col0 to
  # braMin = 1
  # atMin = 1
  # length = 6
  #
  # cur_BrVec <- BrVec[seq(braMin,braMin+length)]
  # cur_AtVec <- AtVec[seq(atMin,atMin+length)]

  #data <- data.frame(at=cur_AtVec, br=cur_BrVec)

  # try sliding arabidopsis from left to right across brassica. Normalise expression by the considered timepoints (in each) mean
  # score inverse of mean distance between points - nb THIS STRETCHES, S.T. POINTS ARE LINED UP, EVEN IF WERE SAMPLED AT DIFFERENT TIMES!!!
  scores <- c()
  num_points <- 10
  for (num_points in seq(min_shift, length(BrVec))) {
    if (num_points <= length(AtVec)) {
      cur_at <- AtVec[seq(length(AtVec) - (num_points-1), length(AtVec))]
      cur_bra <- BrVec[seq(1, num_points)]
    } else {
      cur_at <- AtVec
      cur_bra <- BrVec[seq(1+num_points-length(AtVec), num_points)]
    }

    # rescale both
    cur_at <- cur_at / mean(cur_at)
    cur_bra <- cur_bra / mean(cur_bra)

    score <- 1 / mean((cur_at - cur_bra)^2)
    scores <- c(scores, score)

    adf <- data.frame(expression=cur_at,accession='at',time=seq(1,min(num_points, length(AtVec))))
    bdf <- data.frame(expression=cur_bra,accession='br',time=seq(1,min(num_points, length(AtVec))))
    df <- rbind(adf, bdf)
    # ggplot(df, aes(x=time, y=expression, color=accession))+
    #   geom_point()+
    #   ggtitle(paste0(num_points, '->', score))
  }
  scores[is.na(scores)] <- 0 # if all 0, gives Na
  # best score
  best_ind <- which(scores==max(scores))
  best_num_points <- seq(min_shift, length(BrVec))[best_ind]

  # scores for all shifts (want to check how robust optimal is for key genes):
  all_scores <- data.frame('shift'=seq(min_shift, length(BrVec)), 'score'=scores/sum(scores))
  all_scores$locus_name <- curr_sym

  out <- list(best_num_points, all_scores)
  #return(best_num_points)
  return(out)
}

make_heatmap_quantile <- function(D, title) {
  D$x.sample <- factor(D$x.sample, levels=unique(sort(D$x.sample)))
  D$y.sample <- factor(D$y.sample, levels=unique(sort(D$y.sample)))

  D$quantile <- D$quantile * 100

  p <- ggplot2::ggplot(D)+
    ggplot2::aes(x=x.sample, y=y.sample, fill=quantile) +
    ggplot2::geom_tile()+
    ggplot2::geom_text(ggplot2::aes(label=round(quantile, digits=0)), color='grey', size=1)+
    # viridis::scale_fill_viridis()+
    ggplot2::scale_fill_gradient2(low='royalblue3', mid='white', high='red3', midpoint=50)+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90,
                                                       size=6),
                   axis.text.y = ggplot2::element_text(size=6),
                   plot.title = ggplot2::element_text(hjust=0.5, size=10),
                   legend.position = 'top',
                   legend.justification="right",
                   legend.margin= ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(0,0,-10,-10),
                   legend.text=ggplot2::element_text(size=4, vjust=-0.5),
                   legend.title = ggplot2::element_text(size=8),
                   legend.key.height = ggplot2::unit(0.2, 'cm'),
    )+
    ggplot2::guides(fill=ggplot2::guide_colorbar(label.position='top'))+
    ggplot2::labs(
      x = "",
      y = "",
      title = title
    )

  return(p)
}

get_jobIds <- function(shuffled.data.dir) {
  files <- list.files(shuffled.data.dir)
  tmp <- data.table::tstrsplit(files, '\\.')
  tmp <- data.table::tstrsplit(tmp[[2]], '_')
  jobIds <- unique(paste0(tmp[[2]], '_', tmp[[3]]))
  jobIds <- jobIds[jobIds != 'NA_NA']
  if (length(jobIds) != 1000) {
    print('didnt find 1000 jobIds')
    stop()
  }
  return(jobIds)
}

get_num_registered_genes <- function(file_path) {
  imputed.mean <- readRDS(file_path)
  is.registered.df <- unique(imputed.mean[, c('locus_name', 'is.registered')])
  num.registered <- sum(is.registered.df$is.registered)
  return(num.registered)
}

#c.th <- 0.7
filter.low.variability <- function(exp, c.th) {
  # return gene ids filtered to only include genes with high correlation ( >=c.th)
  # between indiv replicate points, and mean df point

  # calculate correlation
  exp[, C:=stats::cor(mean_cpm, norm.cpm, method='pearson'), by=.(locus_name, accession, tissue)]

  #tmp <- exp[exp$locus_name.model=='MSTRG.8543',]

  # filter to keep if correlationin both accessions > c.th
  exp[, keep:=all(C >= c.th), by=.(locus_name, tissue)]
  keep.df <- exp[exp$keep==TRUE, ]

  #tmp <- keep.df[keep.df$locus_name=='MSTRG.8543',]

  return(unique(keep.df$locus_name))
}


make_heatmap_w_shuffled <- function(D, title) {
  D$x.sample <- factor(D$x.sample, levels=unique(sort(D$x.sample)))
  D$y.sample <- factor(D$y.sample, levels=unique(sort(D$y.sample)))

  p <- ggplot2::ggplot(D)+
    ggplot2::aes(x=x.sample, y=y.sample, fill=log(distance)) +
    ggplot2::geom_tile()+
    #viridis::scale_fill_viridis()+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90,
                                                       size=6),
                   axis.text.y = ggplot2::element_text(size=6),
                   plot.title = ggplot2::element_text(hjust=0.5, size=10),
                   legend.position = 'top',
                   legend.justification="right",
                   legend.margin=ggplot2::margin(0,0,0,0),
                   legend.box.margin = ggplot2::margin(0,0,-10,-10),
                   legend.text=ggplot2::element_text(size=4, vjust=-0.5),
                   legend.title = ggplot2::element_text(size=8),
                   legend.key.height = ggplot2::unit(0.2, 'cm'),
    )+
    ggplot2::guides(fill=ggplot2::guide_colorbar(label.position='top'))+
    viridis::scale_fill_viridis()+
    ggplot2::facet_wrap(~is.shuffled, nrow=1)+
    ggplot2::labs(
      x = "",
      y = "",
      title = title
    )

  return(p)
}

load_shuffled_data <- function(shuffled.data.dir, file.type) {
  # file.type is type of file loading:
  # "comparison": for model.comparison files
  # "mean.sc": for mean_df.sc files
  # "mean_df" : for mean_df files
  # "imputed.mean_df" for them
  # "shifts" for all_shifts

  allowed.types <- c('comparison', 'mean.sc', 'mean_df', 'imputed.mean_df', 'shifts')
  if (!(file.type %in% allowed.types)) {
    print(paste0('file.type must be one of : ', paste0(allowed.types, collapse=', ')))
    stop()
  }
  # format file.type to searchable pattern - handle "."s
  if (file.type=='mean.sc') {
    file.type <- 'mean\\.df\\.sc'
  } else if (file.type=='mean_df') {
    file.type <-  '^mean\\.df'
  } else if (file.type=='imputed.mean_df') {
    file.type <-  'imputed\\.mean\\.df'
  }


  files <- list.files(shuffled.data.dir)
  # make ids
  ids <- get_job_suffixes(shuffled.data.dir)
  i <- 1
  outlist <- rep(list(0), length(ids))
  for (i in 1:length(ids)) {
    id <- ids[i]

    # load files for this job id
    #shift.file <- files[grepl(pattern=paste0('shifts_',id, '\\.'), files)]
    file.idx <- grepl(pattern=paste0(file.type, '_', id, '\\.'), files)
    model.comparison.file <- files[file.idx]

    if (sum(file.idx) != 1) {
      print(paste0('wrong number of files found for : ', file.type, '_', id, '\\.'))
      print(paste0('hits to : ', paste0(model.comparison.file, collapse=' & ')))
      stop()
    }


    #all_shifts <- readRDS(paste0(shuffled.data.dir, shift.file))

    # model comparison only has the best registration used for each gene, so can use
    # to get best shift
    model.comparison <- readRDS(paste0(shuffled.data.dir, model.comparison.file))
    model.comparison$job <- id
    outlist[[i]] <- model.comparison
  }
  out.df <- do.call('rbind', outlist)
  return(out.df)
}

#data.dir <- shuffled.data.dir
get_job_suffixes <- function(data.dir) {
  files <- list.files(data.dir)
  if (length(files)==1) {
    print('error in get_job_suffixes: ')
    print(paste0('Incorrect directory specified : ',  data.dir))
    stop()
  }

  ids <- data.table::tstrsplit(files, '\\.')[[2]]
  ids <- data.table::tstrsplit(ids, '_')
  ids <- unique(paste0(ids[[2]], '_', ids[[3]]))
  ids <- ids[ids != 'NA_NA']

  return(ids)
}


#get_shifted_expression(shift_results, exp)
get_shifted_expression <- function(shift_results, exp) {
  cur_gene <- 'BRAA01G010430.3C'
  shifted_exp <- list()
  for (cur_gene in unique(shift_results$symbol)) {

    # cut to get a single symbol
    test <- exp[exp$locus_name==cur_gene, ]

    # ggplot2::ggplot(test)+
    #   ggplot2::aes(x=timepoint, y=norm.cpm, color=accession)+
    #   ggplot2::geom_point()

    # get the expression vectors for the current gene
    atdf <- test[test$accession=='Col0']
    brdf <- test[test$accession=='Ro18']

    # get the correct times for both
    shift <- shift_results$num.points[shift_results$symbol==cur_gene]
    atdf$timepoint
    atdf$shifted_time <- atdf$timepoint - 7
    atdf$shifted_time <- atdf$shifted_time * 2
    atdf$shifted_time <- atdf$shifted_time + 7
    atdf$shifted_time
    brdf$shifted_time <- brdf$timepoint + (16 - 2*shift)
    brdf$shifted_time
    df <- rbind(atdf, brdf)
    df$shift <- shift

    # normalise each by its mean values in the overlapped time:
    overlapped.times <- unique(intersect(atdf$shifted_time, brdf$shifted_time))
    atdf$norm.cpm <- atdf$norm.cpm / mean(atdf$norm.cpm[atdf$shifted_time %in% overlapped.times])
    brdf$norm.cpm <- brdf$norm.cpm / mean(brdf$norm.cpm[brdf$shifted_time %in% overlapped.times])


    df <- rbind(atdf, brdf)
    df$sc.norm.cpm <- df$norm.cpm
    df$shift <- shift
    # ggplot2::ggplot(df)+
    #   ggplot2::aes(x=shifted_time, y=norm.cpm, color=accession)+
    #   ggplot2::geom_point()

    #adf <- data.frame(accession='Col0', timepoint=at_times, sc.norm.cpm=AtVec, symbol=cur_gene, shift=num_points/3)
    #bdf <- data.frame(accession='Ro18', timepoint=br_times, sc.norm.cpm=BrVec, symbol=cur_gene, shift=num_points/3)
    #df <- rbind(adf, bdf)
    shifted_exp <- c(shifted_exp, list(df))
  }
  shifted_exp <- do.call('rbind', shifted_exp)

  return(shifted_exp)
}


# Commented functions ----

# ro18_rds_file <- '../final_data/rds/ro18_leaf_reannotated.rds'
# sari14_rds_file <- '../final_data/rds/sari14_chiifu_leaf_Jul2020.rds'
#' load_mean_df_sari_ro18 <- function(ro18_rds_file, sari_rds_file) {
#'
#'
#'   ro18.exp <- readRDS(ro18_rds_file)
#'   #ro18.exp$FPKM <- NULL
#'   sari14.exp <- readRDS(sari14_rds_file)
#'   #sari14.exp$TPM <- NULL
#'
#'   # check that comparing the same tissue. Can't easily check
#'   # that the same alignment though.
#'   if (unique(ro18.exp$tissue != unique(sari14.exp$tissue))) {
#'     print("compared tissues aren't the same in sari and exp")
#'     stop()
#'   }
#'
#'   # get the mean of each gene at each timepoint
#'   ro18.exp[, mean_cpm:=mean(norm.cpm), by=.(CDS.model, accession, tissue, timepoint)]
#'   sari14.exp[, mean_cpm:=mean(norm.cpm), by=.(CDS.model, accession, tissue, timepoint)]
#'
#'   exp <- rbind(ro18.exp, sari14.exp)
#'   mean_df <- unique(exp[, c('CDS.model', 'accession', 'tissue', 'timepoint', 'mean_cpm')])
#'   names(mean_df) <-  c('locus_name', 'accession', 'tissue', 'timepoint', 'mean_cpm')
#'   names(exp)[names(exp)=='CDS.model'] <- 'locus_name'
#'   # checking not cutting out key genes - need to know what they're called now
#'   # FT
#'   # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.8543'),]
#'   # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.41993'),]
#'   # # SOC1
#'   # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.15712'),]
#'   # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.26487'),]
#'
#'   # ggplot(tmp, aes(x=timepoint, y=mean_cpm))+
#'   #   geom_line()
#'
#'   # filter mean.df to remove genes with very low expression - remove if max is less than 5, and less than half timepoints expressed
#'   # greater than 1 in both accessions
#'   ro18.df <- mean.df[mean.df$accession=='Ro18']
#'   ro18.df[, keep:=(max(mean_cpm) > 2 | mean(mean_cpm > 1) > 0.5) , by=.(locus_name)]
#'   ro18.keep.genes <- unique(ro18.df$locus_name[ro18.df$keep==TRUE])
#'   #'MSTRG.8543' %in% ro18.keep.genes
#'
#'   sari14.df <- mean.df[mean.df$accession=='sarisha14']
#'   sari14.df[, keep:=(max(mean_cpm) > 2 | mean(mean_cpm > 1) > 0.5) , by=.(locus_name)]
#'   sari14.keep.genes <- unique(sari14.df$locus_name[sari14.df$keep==TRUE])
#'   #'MSTRG.8543' %in% sari14.keep.genes
#'
#'   keep.genes <- intersect(ro18.keep.genes, sari14.keep.genes)
#'   mean.df <- mean.df[mean.df$locus_name %in% keep.genes,]
#'
#'
#'   # filter to remove genes with low correlation between individual timepoints
#'   # and mean timepoints in either accession
#'   keep.genes <- filter.low.variability(exp, 0.7)
#'   #'MSTRG.8543' %in% keep.genes
#'
#'   mean.df <- mean.df[mean.df$locus_name %in% keep.genes, ]
#'   #'MSTRG.8543' %in% mean.df$locus_name
#'
#'   # # checking not cutting out key genes - need to know what they're called now
#'   # # FT
#'   # tmp <- mean.df[mean.df$locus_name %in% c('MSTRG.8543'),]
#'   # tmp <- mean.df[mean.df$locus_name %in% c('MSTRG.41993'),]
#'   # # SOC1
#'   # tmp <- mean.df[mean.df$locus_name %in% c('MSTRG.15712'),]
#'   # tmp <- mean.df[mean.df$locus_name %in% c('MSTRG.26487'),]
#'   #
#'
#'   print(paste0(length(unique(mean.df$locus_name)), ' genes considered in the comparison'))
#'   rm(ro18.df, keep.genes)
#'
#'   exp <- exp[exp$locus_name %in% unique(mean.df$locus_name)]
#'   exp <- subset(exp, select=c('locus_name', 'accession', 'tissue', 'timepoint',
#'                               'norm.cpm', 'group'))
#'   names(exp)[names(exp)=='norm.cpm'] <- 'mean_cpm'
#'   return(list(mean.df, exp))
#' }

# my_scale <- function(V) {
#   V <- (V - mean(V)) / sd(V)
#   return(V)
# }
