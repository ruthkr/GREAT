
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


#file.type <- 'mean_df'
load_shuffled_data <- function(shuffled.data.dir, file.type) {
  # file.type is type of file loading:
  # "comparison": for model.comparison files
  # "mean.sc": for mean_df.sc files
  # "mean_df" : for mean_df files
  # "imputed.mean_df" for them
  # "shifts" for all.shifts

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


    #all.shifts <- readRDS(paste0(shuffled.data.dir, shift.file))

    # model comparison only has the best registration used for each gene, so can use
    # to get best shift
    model.comparison <- readRDS(paste0(shuffled.data.dir, model.comparison.file))
    model.comparison$job <- id
    outlist[[i]] <- model.comparison
  }
  out.df <- do.call('rbind', outlist)
  return(out.df)
}

get_jobIds <- function(shuffled.data.dir) {
  files <- list.files(shuffled.data.dir)
  tmp <- tstrsplit(files, '\\.')
  tmp <- tstrsplit(tmp[[2]], '_')
  jobIds <- unique(paste0(tmp[[2]], '_', tmp[[3]]))
  jobIds <- jobIds[jobIds != 'NA_NA']
  if (length(jobIds) != 1000) {
    print('didnt find 1000 jobIds')
    stop()
  }
  return(jobIds)
}

#data.dir <- shuffled.data.dir
get_job_suffixes <- function(data.dir) {
  files <- list.files(data.dir)
  if (length(files)==1) {
    print('error in get_job_suffixes: ')
    print(paste0('Incorrect directory specified : ',  data.dir))
    stop()
  }

  ids <- tstrsplit(files, '\\.')[[2]]
  ids <- tstrsplit(ids, '_')
  ids <- unique(paste0(ids[[2]], '_', ids[[3]]))
  ids <- ids[ids != 'NA_NA']

  return(ids)
}

get_num_registered_genes <- function(file_path) {
  imputed.mean <- readRDS(file_path)
  is.registered.df <- unique(imputed.mean[, c('locus_name', 'is.registered')])
  num.registered <- sum(is.registered.df$is.registered)
  return(num.registered)
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

  summed.GoIs.df[, scaled.cpm:=my.scale(mean.cpm), by=.(Ara.id, accession)]
  morphology.equiv.df <- data.frame('accession'=c('Col-0', 'R-O-18'), 'floral.transition.time'=c(14, 17))
  summed.GoIs.df$accession <- as.character(summed.GoIs.df$accession)
  summed.GoIs.df$accession[summed.GoIs.df$accession=='Col0'] <- 'Col-0'
  summed.GoIs.df$accession[summed.GoIs.df$accession=='Ro18'] <- 'R-O-18'

  curr.acc <- 'Col-0'
  plot.list <- list()
  for (curr.acc in c('Col-0', 'R-O-18')) {

    curr.p <- ggplot(summed.GoIs.df[summed.GoIs.df$accession==curr.acc,],
                     aes(x=timepoint, y=mean.cpm, color=Ara.name, fill=Ara.name))+
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
  registered.plot.df$accession[registered.plot.df$accession=='Ro18'] <- 'R-O-18'

  p.registered <- ggplot(registered.plot.df, aes(x=shifted_time, y=mean.cpm, color=accession, fill=accession))+
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




#all.rep.data <- all.data.df
scale_all_rep_data <- function(mean_df, all.rep.data, scale.func) {
  # apply the same scaling which done to the mean expression data
  # to all the reps.
  # (subtract mean, and divide by sd), using the values for the mean data
  # as this is what was used to find the best shift.

  # calculate the summary stats to use for the rescaling
  gene.expression.stats <- unique(mean_df[,
                                          .(mean_val=mean(mean.cpm),
                                            sd_val=sd(mean.cpm)),
                                          by=.(locus_name, accession)])

  all.rep.data <- merge(all.rep.data, gene.expression.stats, by=c('locus_name', 'accession'))
  if (scale.func == 'scale') {
    all.rep.data$scaled.norm.cpm <- (all.rep.data$mean.cpm - all.rep.data$mean_val) / all.rep.data$sd_val
  } else if (scale.func == 'my.scale') {
    all.rep.data$scaled.norm.cpm <- (all.rep.data$mean.cpm / all.rep.data$mean_val)
  } else {
    print('invalid scale option for scale_all_rep_data')
    stop()
  }

  out <- subset(all.rep.data, select=c('locus_name', 'accession', 'tissue', 'timepoint',
                                       'scaled.norm.cpm'))

  names(out)[names(out)=='scaled.norm.cpm'] <- 'mean.cpm'

  # ggplot(mean_df[mean_df$locus_name=='BRAA01G000040.3C', ], aes(x=timepoint, y=mean.cpm, color=accession)+
  #          geom_point()


  return(out)
}

my.scale <- function(v) {
  return(v / max(v))
}

make_quantile_heatmaps <- function(Q.mean, Q.scaled, Q.registered) {

  p.mean <- make_heatmap_quantile(Q.mean, 'mean expression')
  p.scaled <- make_heatmap_quantile(Q.scaled, 'scaled mean expression')
  p.registered <- make_heatmap_quantile(Q.registered, 'registered & scaled mean expression')

  p.all <- plot_grid(p.mean, p.scaled, p.registered, ncol=1)

  return(p.all)
}

make_heatmap_quantile <- function(D, title) {
  D$x.sample <- factor(D$x.sample, levels=unique(sort(D$x.sample)))
  D$y.sample <- factor(D$y.sample, levels=unique(sort(D$y.sample)))

  D$quantile <- D$quantile * 100

  p <- ggplot(D, aes(x=x.sample, y=y.sample, fill=quantile))+
    geom_tile()+
    geom_text(aes(label=round(quantile, digits=0)), color='grey', size=1)+
    #scale_fill_viridis()+
    scale_fill_gradient2(low='royalblue3', mid='white', high='red3', midpoint=50)+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90,
                                     size=6),
          axis.text.y = element_text(size=6),
          plot.title = element_text(hjust=0.5, size=10),
          legend.position = 'top',
          legend.justification="right",
          legend.margin=margin(0,0,0,0),
          legend.box.margin = margin(0,0,-10,-10),
          legend.text=element_text(size=4, vjust=-0.5),
          legend.title = element_text(size=8),
          legend.key.height = unit(0.2, 'cm'),
    )+
    guides(fill=guide_colorbar(label.position='top'))+
    ylab('')+
    xlab('')+
    ggtitle(title)

  return(p)
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

ecdf_fun <- function(x,perc) {
  # estimate the quantile of the perc value in the x distribution
  # https://stats.stackexchange.com/questions/50080/estimate-quantile-of-value-in-a-vector/114493
  ecdf(x)(perc)
}


make_shuffled_data_heatmaps <- function(D.mean, D.scaled, D.registered) {

  p.mean <- make_heatmap_w_shuffled(D.mean, 'mean expression')
  p.scaled <- make_heatmap_w_shuffled(D.scaled, 'scaled mean expression')
  p.registered <- make_heatmap_w_shuffled(D.registered, 'registered & scaled mean expression')

  p.all <- plot_grid(p.mean, p.scaled, p.registered, ncol=1)

  return(p.all)
}

make_heatmap_w_shuffled <- function(D, title) {
  D$x.sample <- factor(D$x.sample, levels=unique(sort(D$x.sample)))
  D$y.sample <- factor(D$y.sample, levels=unique(sort(D$y.sample)))

  p <- ggplot(D, aes(x=x.sample, y=y.sample, fill=log(distance)))+
    geom_tile()+
    #scale_fill_viridis()+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90,
                                     size=6),
          axis.text.y = element_text(size=6),
          plot.title = element_text(hjust=0.5, size=10),
          legend.position = 'top',
          legend.justification="right",
          legend.margin=margin(0,0,0,0),
          legend.box.margin = margin(0,0,-10,-10),
          legend.text=element_text(size=4, vjust=-0.5),
          legend.title = element_text(size=8),
          legend.key.height = unit(0.2, 'cm'),
    )+
    guides(fill=guide_colorbar(label.position='top'))+
    scale_fill_viridis()+
    facet_wrap(~is.shuffled, nrow=1)+
    ylab('')+
    xlab('')+
    ggtitle(title)

  return(p)
}

load_shuffled_distances <- function(dir.name, D.type, numJobs, numRuns) {
  # dir.name the directory to look in (e.g. "rescaled_shuffled_distance") in intermediate_data/gene_registration
  # D.type in c('D.mean', 'D.scaled', 'D.registered') depending
  # how the shuffled expression data was transformed
  # mean of biological replicates
  # that + expression scaled(center=T, scale=T)
  # that + registered (scale using compared points only)

  # numJobs : how many seperate jobs were submitted
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
  ro18.exp[, mean.cpm:=mean(norm.cpm), by=.(CDS.model, accession, tissue, timepoint)]
  sari14.exp[, mean.cpm:=mean(norm.cpm), by=.(CDS.model, accession, tissue, timepoint)]

  exp <- rbind(ro18.exp, sari14.exp)
  mean_df <- unique(exp[, c('CDS.model', 'accession', 'tissue', 'timepoint', 'mean.cpm')])
  names(mean_df) <-  c('locus_name', 'accession', 'tissue', 'timepoint', 'mean.cpm')
  names(exp)[names(exp)=='CDS.model'] <- 'locus_name'
  # checking not cutting out key genes - need to know what they're called now
  # FT
  # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.8543'),]
  # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.41993'),]
  # # SOC1
  # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.15712'),]
  # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.26487'),]

  # ggplot(tmp, aes(x=timepoint, y=mean.cpm))+
  #   geom_line()

  # filter mean_df to remove genes with very low expression - remove if max is less than 5, and less than half timepoints expressed
  # greater than 1 in both accessions
  ro18.df <- mean_df[mean_df$accession=='Ro18']
  ro18.df[, keep:=(max(mean.cpm) > 2 | mean(mean.cpm > 1) > 0.5) , by=.(locus_name)]
  ro18.keep.genes <- unique(ro18.df$locus_name[ro18.df$keep==TRUE])
  #'MSTRG.8543' %in% ro18.keep.genes

  sari14.df <- mean_df[mean_df$accession=='sarisha14']
  sari14.df[, keep:=(max(mean.cpm) > 2 | mean(mean.cpm > 1) > 0.5) , by=.(locus_name)]
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
  names(exp)[names(exp)=='norm.cpm'] <- 'mean.cpm'
  return(list(mean_df, exp))
}

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
#'   ro18.exp[, mean.cpm:=mean(norm.cpm), by=.(CDS.model, accession, tissue, timepoint)]
#'   sari14.exp[, mean.cpm:=mean(norm.cpm), by=.(CDS.model, accession, tissue, timepoint)]
#'
#'   exp <- rbind(ro18.exp, sari14.exp)
#'   mean_df <- unique(exp[, c('CDS.model', 'accession', 'tissue', 'timepoint', 'mean.cpm')])
#'   names(mean_df) <-  c('locus_name', 'accession', 'tissue', 'timepoint', 'mean.cpm')
#'   names(exp)[names(exp)=='CDS.model'] <- 'locus_name'
#'   # checking not cutting out key genes - need to know what they're called now
#'   # FT
#'   # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.8543'),]
#'   # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.41993'),]
#'   # # SOC1
#'   # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.15712'),]
#'   # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.26487'),]
#'
#'   # ggplot(tmp, aes(x=timepoint, y=mean.cpm))+
#'   #   geom_line()
#'
#'   # filter mean_df to remove genes with very low expression - remove if max is less than 5, and less than half timepoints expressed
#'   # greater than 1 in both accessions
#'   ro18.df <- mean_df[mean_df$accession=='Ro18']
#'   ro18.df[, keep:=(max(mean.cpm) > 2 | mean(mean.cpm > 1) > 0.5) , by=.(locus_name)]
#'   ro18.keep.genes <- unique(ro18.df$locus_name[ro18.df$keep==TRUE])
#'   #'MSTRG.8543' %in% ro18.keep.genes
#'
#'   sari14.df <- mean_df[mean_df$accession=='sarisha14']
#'   sari14.df[, keep:=(max(mean.cpm) > 2 | mean(mean.cpm > 1) > 0.5) , by=.(locus_name)]
#'   sari14.keep.genes <- unique(sari14.df$locus_name[sari14.df$keep==TRUE])
#'   #'MSTRG.8543' %in% sari14.keep.genes
#'
#'   keep.genes <- intersect(ro18.keep.genes, sari14.keep.genes)
#'   mean_df <- mean_df[mean_df$locus_name %in% keep.genes,]
#'
#'
#'   # filter to remove genes with low correlation between individual timepoints
#'   # and mean timepoints in either accession
#'   keep.genes <- filter.low.variability(exp, 0.7)
#'   #'MSTRG.8543' %in% keep.genes
#'
#'   mean_df <- mean_df[mean_df$locus_name %in% keep.genes, ]
#'   #'MSTRG.8543' %in% mean_df$locus_name
#'
#'   # # checking not cutting out key genes - need to know what they're called now
#'   # # FT
#'   # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.8543'),]
#'   # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.41993'),]
#'   # # SOC1
#'   # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.15712'),]
#'   # tmp <- mean_df[mean_df$locus_name %in% c('MSTRG.26487'),]
#'   #
#'
#'   print(paste0(length(unique(mean_df$locus_name)), ' genes considered in the comparison'))
#'   rm(ro18.df, keep.genes)
#'
#'   exp <- exp[exp$locus_name %in% unique(mean_df$locus_name)]
#'   exp <- subset(exp, select=c('locus_name', 'accession', 'tissue', 'timepoint',
#'                               'norm.cpm', 'group'))
#'   names(exp)[names(exp)=='norm.cpm'] <- 'mean.cpm'
#'   return(list(mean_df, exp))
#' }

#c.th <- 0.7
filter.low.variability <- function(exp, c.th) {
  # return gene ids filtered to only include genes with high correlation ( >=c.th)
  # between indiv replicate points, and mean df point

  # calculate correlation
  exp[, C:=cor(mean.cpm, norm.cpm, method='pearson'), by=.(locus_name, accession, tissue)]

  #tmp <- exp[exp$locus_name.model=='MSTRG.8543',]

  # filter to keep if correlationin both accessions > c.th
  exp[, keep:=all(C >= c.th), by=.(locus_name, tissue)]
  keep.df <- exp[exp$keep==TRUE, ]

  #tmp <- keep.df[keep.df$locus_name=='MSTRG.8543',]

  return(unique(keep.df$locus_name))
}

load_mean_df <- function() {

  #setwd('/Volumes/Research-Projects/bravo/alex/BRAVO_rna-seq/scripts/')
  rds_file <- 'ro18_chiifu_apex' # don't include the .rds # the name of the brassica data to load
  sumBrassicas <- F # if false use seperate brassica genes, and compare to repeated arabidopsis genes. If true, sume copies of each brassica and compare to arabidopsis
  #datapath <- paste0('../final_data/rds/', rds_file, '.rds')


  #### specify the genes to be used in the comparison: ####
  # GoIs <- fread(paste0('../graphs/', rds_file, '/comparison_genes.tsv')) # read the list of Arabidopsis id's for the genes of interest

  # ------ changed by Ruth
  GoIs <- fread(paste0('graphs/', rds_file, '/comparison_genes.tsv'))

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
  mean_df <- unique(exp[, c('locus_name', 'accession', 'tissue', 'timepoint', 'mean.cpm')])



  # filter mean_df to remove genes with very low expression - remove if max is less than 5, and less than half timepoints expressed
  # greater than 1
  ro18.df <- mean_df[mean_df$accession=='Ro18']
  ro18.df[, keep:=(max(mean.cpm) > 2 | mean(mean.cpm > 1) > 0.5) , by=.(locus_name)]
  keep.genes <- unique(ro18.df$locus_name[ro18.df$keep==TRUE])
  mean_df <- mean_df[mean_df$locus_name %in% keep.genes,]
  print(paste0(length(unique(mean_df$locus_name)), ' genes considered in the comparison'))
  rm(ro18.df, keep.genes)

  exp <- exp[exp$locus_name %in% unique(mean_df$locus_name)]
  exp <- subset(exp, select=c('locus_name', 'accession', 'tissue', 'timepoint',
                              'norm.cpm', 'group'))
  names(exp)[names(exp)=='norm.cpm'] <- 'mean.cpm'
  return(list(mean_df, exp))
}

# calculate the comparison stats
calc.AIC <- function(logL, num.params) {
  return((-2*logL) + 2*num.params)
}

calc.BIC <- function(logL, num.params, num.obs) {
  return((-2*logL) + log(num.obs) * num.params)
}


change.accession.names <- function(mean_df, all.data.df, transformed.timecourse) {
  # set the "transformed.timecourse" accession to "Col0", and the other one to "Ro18"

  # error checking
  if (length(unique(mean_df$accession)) != 2) {
    stop('Error in change.accession.names() : comparison must be made between two accessions!')
  }

  # store these, to rename at the end
  original.transformed.timecourse.name <- transformed.timecourse
  original.other.accession.name <- as.character(unique(mean_df$accession[mean_df$accession!=transformed.timecourse]))

  # change mean_df
  new.mean_df.accession <- mean_df$accession
  new.mean_df.accession[mean_df$accession==transformed.timecourse] <- 'Col0'
  new.mean_df.accession[mean_df$accession!=transformed.timecourse] <- 'Ro18'
  mean_df$accession <- new.mean_df.accession

  # change all.data.df
  new.all.data.df.accession <- all.data.df$accession
  new.all.data.df.accession[all.data.df$accession==transformed.timecourse] <- 'Col0'
  new.all.data.df.accession[all.data.df$accession!=transformed.timecourse] <- 'Ro18'
  all.data.df$accession <- new.all.data.df.accession

  return(list('mean_df'=mean_df,
              'all.data.df'=all.data.df,
              'original.transformed.accession.name'=original.transformed.timecourse.name,
              'original.other.accession.name'=original.other.accession.name))
}


# stretches=c(1, 1.5, 2.0)
# initial.rescale=FALSE
# do_rescale=FALSE
# min_num_overlapping_points = 4
# test.genes <- unique(mean_df$locus_name)[1:101]
# #test.genes <- 'MSTRG.10244'
# mean_df <- mean_df[mean_df$locus_name %in% test.genes,]
# all.data.df <- all.data.df[all.data.df$locus_name %in% test.genes,]
# shift_extreme=4
# transformed.timecourse <- 'Ro18'

prepare_scaled_and_registered_data <- function(mean_df, all.data.df, stretches,
                                               initial.rescale, do_rescale,
                                               min_num_overlapping_points, shift_extreme,
                                               transformed.timecourse) {
  ## APPLY NORMALISATION OF EXPRESSION FOR EACH GENE ACROSS ALL TIMEPOINTS ##

  # hardcoded all the functions to use 'Col0', and 'Ro18'. Rather than fix all instances
  # of that, just temporarily rename them here, and turn back at the end.
  # will apply stretch, and shift to the "transformed.timecourse" accession (which should be the quicker one...)
  L <- change.accession.names(mean_df, all.data.df, transformed.timecourse)
  mean_df <- L[['mean_df']]
  all.data.df <- L[['all.data.df']]
  original.transformed.accession <- L[['original.transformed.accession.name']]
  original.other.accession <- L[['original.other.accession.name']]



  mean_df.sc <- copy(mean_df)
  # specify what kind of scaling
  mean_df.sc[, sc.mean.cpm:=scale(mean.cpm, scale=TRUE, center=TRUE), by=.(locus_name, accession)]
  #mean_df.sc[, sc.mean.cpm:=my.scale(mean.cpm), by=.(locus_name, accession)]

  ## APPLY INDIVIDUAL SHIFTING ######
  # optimise transformations applied to arabidopsis gene profile to map onto the brassicas - shift in x direction, using mean for mapping,
  # and only rescale expression using the candidate timepoints in common.

  # specify which data to use for registering.
  # whether prior rescaled mean, or mean data should be used for registration
  if (initial.rescale==TRUE) {
    # apply rescale to mean_df prior to registration
    to.shift.df <- copy(mean_df.sc)
    to.shift.df$mean.cpm <- to.shift.df$sc.mean.cpm
    to.shift.df$sc.mean.cpm <- NULL

    # apply THE SAME rescale to all.data.df prior to registration
    #all.data.df <- scale_all_rep_data(mean_df, all.data.df, 'my.scale')
    all.data.df <- scale_all_rep_data(mean_df, all.data.df, 'scale')


    # sanity plot that rescale all data worked
    # ggplot(all.data.df[all.data.df$locus_name=='BRAA01G000040.3C'], aes(x=timepoint, y=mean.cpm, color=accession))+
    #   geom_point()
  } else {
    to.shift.df <- copy(mean_df)
  }
  # ggplot(to.shift.df[to.shift.df$locus_name=='BRAA03G004600.3C'], aes(x=timepoint, y=mean.cpm, color=accession))+
  #   geom_point()
  # tst <- all.data.df
  # ggplot(tst[tst$locus_name=='BRAA01G000040.3C'], aes(x=timepoint, y=mean.cpm, color=accession))+
  #   geom_point()

  # calculate the best registration. Returns all tried registrations, best stretch and shift combo,
  # and AIC/BIC stats for comparison of best registration model to seperate models for expression of
  # each gene in Ro18 and Col0.
  L <- get_best_stretch_and_shift(to.shift.df, all.data.df, stretches, do_rescale, min_num_overlapping_points, shift_extreme)
  all_shifts <- L[['all_shifts']]
  best_shifts <- L[['best_shifts']]
  model.comparison.dt <- L[['model.comparison.dt']]


  # report model comparison results
  model.comparison.dt$BIC.registered.is.better <- (model.comparison.dt$registered.BIC < model.comparison.dt$seperate.BIC)
  model.comparison.dt$AIC.registered.is.better <- (model.comparison.dt$registered.AIC < model.comparison.dt$seperate.AIC)
  model.comparison.dt$ABIC.registered.is.better <- (model.comparison.dt$BIC.registered.is.better & model.comparison.dt$AIC.registered.is.better)
  print('################## Model comparison results #######################')
  print(paste0('AIC finds registration better than seperate for :', sum(model.comparison.dt$AIC.registered.is.better), ' / ', nrow(model.comparison.dt)))
  print(paste0('BIC finds registration better than seperate for :', sum(model.comparison.dt$BIC.registered.is.better), ' / ', nrow(model.comparison.dt)))
  print(paste0('AIC & BIC finds registration better than seperate for :', sum(model.comparison.dt$ABIC.registered.is.better), ' / ', nrow(model.comparison.dt)))
  print('###################################################################')


  # get the best-shifted and stretched mean gene expression, only to genes which registration is better than
  # seperate models by BIC. Don't stretch out, or shift genes for which seperate is better.
  # registration is applied to col0.
  shifted.mean_df <- apply_shift_to_registered_genes_only(to.shift.df, best_shifts, model.comparison.dt)
  # shifted.mean_df <- apply_best_shift(to.shift.df, best_shifts) # can be NA if exactly tied for what the best shift was

  # GOI <- 'MSTRG.11237'
  # ggplot(all.data.df[all.data.df$locus_name==GOI],
  #        aes(x=timepoint, y= mean.cpm, color=accession))+
  #   geom_point()
  #
  # #sanity plot that done right
  # ggplot(shifted.mean_df[shifted.mean_df$locus_name==GOI],
  #        aes(x=shifted_time, y=mean.cpm, color=accession))+
  #   geom_point()+
  #   geom_line()


  # impute arabidopsis values at times == to the observed brassica points for each shifted arabidopsis gene
  # so can compare using heatmap.
  # arabidopsis curves are the ones that been shifted around. Linear impute values for these
  # curves so that brassica samples can be compared to an arabidopsis point.
  imputed.mean_df <- impute_arabidopsis_values(shifted.mean_df)

  #sanity plot that done right
  # ggplot(shifted.mean_df[shifted.mean_df$locus_name=='BRAA01G001540.3C'],
  #        aes(x=shifted_time, y=mean.cpm, color=accession))+
  #   geom_point()+
  #   geom_line()
  # ggplot(imputed.mean_df[imputed.mean_df$locus_name=='BRAA01G001540.3C'],
  #        aes(x=shifted_time, y=mean.cpm, color=accession))+
  #   geom_point()+
  #   geom_line()


  # fix the accession names to the ones actually passed in:
  mean_df <- fix.accessions(mean_df, original.transformed.accession, original.other.accession)
  mean_df.sc <- fix.accessions(mean_df.sc, original.transformed.accession, original.other.accession)
  imputed.mean_df <- fix.accessions(imputed.mean_df, original.transformed.accession, original.other.accession)


  OUT <- list('mean_df'=mean_df,
              'mean_df.sc'=mean_df.sc,
              'imputed.mean_df'=imputed.mean_df,
              'all.shifts'=all_shifts,
              'model.comparison'=model.comparison.dt)
}

fix.accessions <- function(df, original.transformed.accession, original.other.accession) {
  # swap Col0 with original.transformed.accession, and Ro18 with original.other.accession
  new.df.accession <- df$accession
  new.df.accession[df$accession=='Col0'] <- original.transformed.accession
  new.df.accession[df$accession=='Ro18'] <- original.other.accession
  df$accession <- new.df.accession

  return(df)
}


get_best_result <- function(df) {
  # return TRUE/FALSE vector. TRUE for the smallest score
  # if tied for this, true for the one with the smallest stretch. (1x is smaller than 0.75x though)
  # if tied, then the one with the smallest shift

  is.best <- df$score==min(df$score)
  if(sum(is.best)==1) {
    return(is.best)
  } else {
    cand.stretches <- df$stretch[is.best]
    # get the stretch with the best score, with the smallest divergence from 1
    min.stretch <- unique(cand.stretches[abs(cand.stretches-1)==min(abs(cand.stretches-1))])
    is.best[df$stretch != min.stretch] <- FALSE
    if(sum(is.best)==1) {
      return(is.best)
    } else {
      cand.shifts <- df$shift[is.best]
      min_shift <- unique(cand.shifts[cand.shifts==min(cand.shifts)])
      is.best[df$shift != min_shift] <- FALSE
      if (sum(is.best)==1) {
        return(is.best)
      } else {
        stop('error in get_best_result, somehow STILL more than one best shift tied..?')
      }
    }
  }
}

# to.shift.df
# all.data.df
# stretches
# do_rescale
# min_num_overlapping_points
# shift_extreme
get_best_stretch_and_shift <- function(to.shift.df, all.data.df, stretches, do_rescale, min_num_overlapping_points, shift_extreme) {
  # for each stretch in stretches, calculates best shift, by comparing SUM of squares difference.
  # for the best shift in each stretch, compares to seperate models to calculate AIC/BIC under registration,
  # / no registration.

  # returns:
  # - all_shifts : all the combos of stretching and shifting tried for each gene
  # - best_shifts : the best stretch and shift combo found for each gene, as well as info for scaling etc
  # - model_comparison.dt : AIC / BIC scores for best registerd model found, compared to seperate model for
  # each genes expression in the 2 accessions.


  if (!('Col0' %in% all.data.df$accession & 'Ro18' %in% all.data.df$accession)) {
    stop('get_best_stretch_and_shift() : all.data.df accessions should have been
         converted to Col0 and Ro18.')
  }

  all_all_shifts <- rep(list(0), length(stretches))
  all_best_shifts <- rep(list(0), length(stretches))
  all_model_comparison.dt <- rep(list(0), length(stretches))
  i <- 3
  for (i in 1:length(stretches)) {
    stretch <- stretches[i]
    print(paste0('testing models for stretch factor = ', stretch))
    # calculate all the shift scores given this stretch. Score is mean(dist^2), over overlapping points
    # if do_rescale=T, is rescaled by the mean FOR THE OVERLAPPING POINTS. (but not by the SD.)

    # ggplot(to.shift.df[to.shift.df$locus_name=='MSTRG.12467',], aes(x=timepoint, y=mean.cpm, color=accession))+
    #   geom_point()

    all_shifts <- calculate_all_best_shifts(to.shift.df, stretch_factor=stretch, do_rescale, min_num_overlapping_points, shift_extreme)

    all_shifts <- unique(all_shifts) # ensure no dulicated rows

    # cut down to single best shift for each gene
    all_shifts[, is.best:=get_best_result(.SD), by=.(gene)]
    best_shifts <- all_shifts[is.best==TRUE,]
    all_shifts$is.best <- NULL

    if (nrow(best_shifts)!= length(unique(all.data.df$locus_name))) {
      stop('get_best_stretch_and_shift() : got non-unique best shifts in best_shifts')
    }

    # calculate the BIC & AIC for the best shifts found with this stretch.compared to treating the
    # gene's expression seperately in Col0 and Ro18
    model.comparison.dt <- calculate_all_model_comparison_stats(all.data.df, best_shifts)
    # add info on the stretch and shift applied
    model.comparison.dt <- merge(model.comparison.dt, best_shifts[, c('gene', 'stretch', 'shift'),],
                                 by='gene')

    # record the results for the current stretch factor
    all_all_shifts[[i]] <- all_shifts
    all_best_shifts[[i]] <- best_shifts
    all_model_comparison.dt[[i]] <- model.comparison.dt
  }

  all_shifts <- do.call('rbind', all_all_shifts) # all the combinations of shift, and stretch tried
  all_best_shifts <- do.call('rbind', all_best_shifts) # the best shifts for each stretch
  all_model_comparison.dt <- do.call('rbind', all_model_comparison.dt) # model comparison of best shift (for each stretch) to seperate modeles

  # get the best registration applied (best stretch, and best shift) for each gene
  # picking by bic alone will favour fewer overlapping (considered) datapoints. Pick best in order to maximise
  # how much better register.BIC is than seperate.BIC
  all_model_comparison.dt$delta.BIC <- all_model_comparison.dt$registered.BIC - all_model_comparison.dt$seperate.BIC
  all_model_comparison.dt[, is.best:=(delta.BIC==min(delta.BIC)), by=.(gene)] # best is one for which registered.BIC is as small as possible compared to seperate.BIC
  best_model_comparison.dt <- all_model_comparison.dt[all_model_comparison.dt$is.best==TRUE]
  # if there's a tie for best registration for a gene, keep the first one as the best
  if (any(duplicated(best_model_comparison.dt$gene))) {
    print(paste0('found ', sum(duplicated(best_model_comparison.dt$gene)), ' tied optimal registrations. Removing dupliates'))
    best_model_comparison.dt <- best_model_comparison.dt[!(duplicated(best_model_comparison.dt$gene)),]
  }
  best_model_comparison.dt$delta.BIC <- NULL
  # cut down best shifts to the best shift for the best stretch only
  best_shifts <- merge(all_best_shifts, best_model_comparison.dt[, c('gene', 'stretch', 'shift')],
                       by=c('gene', 'stretch', 'shift'))
  stopifnot(nrow(best_shifts)==length(unique(to.shift.df$locus_name))) # there should be only 1 best shift for each gene

  return(list('all_shifts'=all_shifts,
              'best_shifts'=best_shifts,
              'model.comparison.dt'=best_model_comparison.dt))
}


apply_shift_to_registered_genes_only <- function(to.shift.df, best_shifts, model.comparison.dt) {

  # genes for which registration model is better than seperate model
  genes.to.register <- model.comparison.dt$gene[model.comparison.dt$BIC.registered.is.better]
  # apply the registration transformation to these genes
  if (length(genes.to.register > 0)) {
    register.dt <- to.shift.df[to.shift.df$locus_name %in% genes.to.register,]
    registered.dt <- apply_best_shift(register.dt, best_shifts)
    registered.dt$is.registered <- TRUE
  }

  # genes for which the seperate model is better than registration model
  genes.to.keep.seperate <- model.comparison.dt$gene[!(model.comparison.dt$BIC.registered.is.better)]
  # generate the columns for these needed to concat. with registered.dt
  seperate.dt <- to.shift.df[to.shift.df$locus_name %in% genes.to.keep.seperate,]
  seperate.dt$stretched.time.delta <- 0 # in order to ensure that seperate copy
  # print('line 594')
  # print(min(timepoint))
  seperate.dt[, stretched.time.delta:=timepoint - min(timepoint), by=.(locus_name, accession)]
  seperate.dt$shifted_time <- seperate.dt$stretched.time.delta + 11 # add eleven, as this is done for the registered genes
                                                                    # to make comparible between Ro18 and Col. Therefore need to to this
                                                                    # here, to keep unregistered col0 in same frame as
                                                                    # stretch 1, shift 0 registered genes.
  seperate.dt$is.registered <- FALSE

  if (length(genes.to.register > 0)) {
    out.dt <- rbind(registered.dt, seperate.dt)
  } else {
    out.dt <- seperate.dt
  }

  return(out.dt)
}

calculate_all_model_comparison_stats <- function(all.data.df, best_shifts) {
  # wrapper to apply compare_registered_to_unregistered_model() to all the genes

  if (!('Col0' %in% unique(all.data.df$accession) &
                           'Ro18' %in% unique(all.data.df$accession))) {
    stop("error in calculate_all_model_comparison_stats() :
         all.data.df doesn't have the correct accession info - should have been
         converted to Ro18 & Col0")
  }


  # apply the registration to the all rep data, so can use for model comparison.


  shifted.all.data.df <- apply_best_shift(all.data.df, best_shifts)

  # # check that have now got scaled, stretched time for both genes.
  # tst <- all.data.df
  # ggplot(tst[tst$locus_name=='BRAA01G000040.3C'], aes(x=shifted_time, y=mean.cpm, color=accession))+
  #   geom_point()
  # ggplot(tst[tst$locus_name=='BRAA01G000040.3C'], aes(x=timepoint, y=mean.cpm, color=accession))+
  #   geom_point()

  print('calculating registration vs different expression comparison AIC & BIC...')

  genes <- unique(shifted.all.data.df$locus_name)
  out.sepAIC <- rep(0, length(genes))
  out.combAIC <- rep(0, length(genes))
  out.sepBIC <- rep(0, length(genes))
  out.combBIC <- rep(0, length(genes))

  i <- 1
  #i <- sample(1:length(genes), 1)
  for(i in 1:length(genes)) {
    if( i %% 100 == 0) {
      print(paste0(i, ' / ', length(genes)))
    }

    #curr.sym <- 'BRAA06G010220.3C'
    curr.sym <- genes[i]
    L <- compare_registered_to_unregistered_model(curr.sym, shifted.all.data.df, is.testing=FALSE)
    out.sepAIC[i] <- L[[1]]
    out.combAIC[i] <- L[[2]]
    out.sepBIC[i] <- L[[3]]
    out.combBIC[i] <- L[[4]]
  }

  out <- data.table(data.frame('gene'=genes, 'seperate.AIC'=out.sepAIC, 'registered.AIC'=out.combAIC,
                               'seperate.BIC'=out.sepBIC, 'registered.BIC'=out.combBIC))
  return(out)
}

# curr.sym <- 'BRAA01G001320.3C'
# all.data.df <- all.data.df
# is.testing <- TRUE
compare_registered_to_unregistered_model <- function(curr.sym, all.data.df, is.testing) {
  # compare the overlapping timepoints in brassica and arabidopsis after the best registration,
  # and without registration (use the same timepoints for both models).
  # use the stretched data for both models, whether considering as registered or not.


  curr.data.df <- all.data.df[all.data.df$locus_name==curr.sym]

  # print('line 662')
  # print(curr.data.df)

  # flag the timepoints to be used in the modelling, only the ones which overlap!
  curr.data.df <- get_compared_timepoints(curr.data.df)

  # ggplot(curr.data.df, aes(x=shifted_time, y=mean.cpm, shape=is.compared, color=accession))+
  #   geom_point()

  # cut down to the data for each model
  ara.spline.data <- curr.data.df[curr.data.df$is.compared==TRUE &
                                    curr.data.df$accession=='Col0', ]
  bra.spline.data <- curr.data.df[curr.data.df$is.compared==TRUE &
                                    curr.data.df$accession=='Ro18', ]
  combined.spline.data <- curr.data.df[curr.data.df$is.compared==TRUE, ]

  # fit the models - fit regression splines.
  # http://www.utstat.utoronto.ca/reid/sta450/feb23.pdf
  # for cubic spline, K+3 params where K=num.knots
  # as can omit constant term
  num.spline.params <- 6 # number of parameters for each spline fitting (degree and this used to calculate num knots).
  num.registration.params <- 2 # stretch, shift
  num.obs <- nrow(combined.spline.data)

  # print('line 676')
  # print(ara.spline.data)
  # print(bra.spline.data)


  ara.fit <- lm(mean.cpm~splines::bs(shifted_time, df=num.spline.params, degree=3), data=ara.spline.data)
  bra.fit <- lm(mean.cpm~splines::bs(shifted_time, df=num.spline.params, degree=3), data=bra.spline.data)
  combined.fit <- lm(mean.cpm~splines::bs(shifted_time, df=num.spline.params, degree=3), data=combined.spline.data)
  # calculate the log likelihoods
  ara.logLik <- logLik(ara.fit)
  bra.logLik <- logLik(bra.fit)
  seperate.logLik <- ara.logLik + bra.logLik # logLikelihoods, so sum
  combined.logLik <- logLik(combined.fit)

  # calculate the comparison.stats - - AIC, BIC, smaller is better!
  # 2*num.spline.params as fitting seperate models for Ara * Col
  seperate.AIC <- calc.AIC(seperate.logLik, 2*num.spline.params)
  combined.AIC <- calc.AIC(combined.logLik, num.spline.params+num.registration.params)

  seperate.BIC <- calc.BIC(seperate.logLik, 2*num.spline.params, num.obs)
  combined.BIC <- calc.BIC(combined.logLik, num.spline.params+num.registration.params, num.obs)


  if (is.testing==TRUE) {
    ara.pred <- predict(ara.fit)
    ara.pred.df <- unique(data.frame('shifted_time'=ara.spline.data$shifted_time,
                                     'mean.cpm'=ara.pred, 'accession'='Col0'))
    bra.pred <- predict(bra.fit)
    bra.pred.df <- unique(data.frame('shifted_time'=bra.spline.data$shifted_time,
                                     'mean.cpm'=bra.pred, 'accession'='Ro18'))

    combined.pred <- predict(combined.fit)
    combined.pred.df <- unique(data.frame('shifted_time'=combined.spline.data$shifted_time,
                                          'mean.cpm'=combined.pred, 'accession'='registered'))
    spline.df <- rbind(ara.pred.df, bra.pred.df, combined.pred.df)

    ggplot(data=combined.spline.data, aes(x=shifted_time, y=mean.cpm,
                                          colour=accession))+
      geom_point()+
      geom_line(data=spline.df)+
      ggtitle(paste0(curr.sym, ' : sep AIC:combo AIC=', round(seperate.AIC), ':', round(combined.AIC),
                     ', sep BIC: combo BIC=', round(seperate.BIC), ':', round(combined.BIC)))
    ggsave(paste0('./testing/fitted_splines/', curr.sym, '_', max(ara.pred.df$shifted_time), '.pdf'))
  }

  return(list(seperate.AIC, combined.AIC, seperate.BIC,combined.BIC))
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



make_heatmap <- function(D, ylabel, y.axis.fontsize=6) {
  D$x.sample <- factor(D$x.sample, levels=unique(sort(D$x.sample)))
  D$y.sample <- factor(D$y.sample, levels=unique(sort(D$y.sample)))

  p <- ggplot(D, aes(x=x.sample, y=y.sample, fill=log(distance)))+
    geom_tile()+
    #scale_fill_viridis()+
    #theme_classic()+
    theme(axis.text.x = element_text(angle=90,
                                     size=6),
          axis.text.y = element_text(size=y.axis.fontsize),
          plot.title = element_text(hjust=0.5, size=10),
          plot.margin = margin(0,0,-10,0),
          panel.background = element_blank(),
          legend.position = 'top',
          legend.justification="right",
          legend.margin=margin(0,0,0,0),
          legend.box.margin = margin(0,0,-10,-10),
          legend.text=element_text(size=4, vjust=-0.5),
          legend.title = element_text(size=8),
          legend.key.height = unit(0.2, 'cm'),
    )+
    guides(fill=guide_colorbar(label.position='top'))+
    scale_fill_viridis()+
    facet_wrap(~title, nrow=1)+
    ylab(ylabel)+
    xlab('')
  #ggtitle(title)

  return(p)
}


# mean_df <- real.mean_df
# mean_df.sc <- real.sc.df
# imputed.mean_df <- imputed.mean_df
calculate_between_sample_distance <- function(mean_df, mean_df.sc, imputed.mean_df) {

  ### convert all to wide format ready for distance calculation

  # mean_df
  sample.id.cols <- c('accession','timepoint')
  gene.col <- c('locus_name')
  expression.col <- 'mean.cpm'
  mean.dt.w <- reformat_for_distance_calculation(mean_df, sample.id.cols, gene.col, expression.col)

  # normalised mean_df
  sample.id.cols <- c('accession','timepoint')
  gene.col <- c('locus_name')
  expression.col <- c('sc.mean.cpm')
  mean.dt.sc.w <- reformat_for_distance_calculation(mean_df.sc, sample.id.cols, gene.col, expression.col)

  # imputed.mean_df - all genes
  sample.id.cols <- c('accession','shifted_time')
  gene.col <- c('locus_name')
  expression.col <- c('mean.cpm')
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


# unrandom.mean_df <- mean_df
# unrandom.all.df <- all.data.df
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

shuffle_ro18_gene_names <- function(mean_df, out.all.df) {
  # shuffle the identities of the genes in the brassica
  # can't just do shuffle, becuase need to preserve which timepoints are from the same gene
  out.mean_df <- copy(mean_df)
  out.all.df <- copy(out.all.df)

  # make the gene lookup table for the same shuffled genes for both
  brassica.genes <- unique(out.mean_df$locus_name[out.mean_df$accession=='Ro18'])
  shuffled.genes <- sample(brassica.genes)
  shuffle.gene.lookup <- data.table(data.frame('gene.id'=brassica.genes, 'shuffled.id'=shuffled.genes))

  # change the gene names for the mean_df
  out.mean_df <- swap_gene_names(out.mean_df, shuffle.gene.lookup)
  # change the gene names for the all.df
  out.all.df <- swap_gene_names(out.all.df, shuffle.gene.lookup)

  return(list(out.mean_df, out.all.df))
}

swap_gene_names <- function(df, shuffle.gene.lookup) {
  replacement.genes <- sapply(df$locus_name[df$accession=='Ro18'],
                              function(x) shuffle.gene.lookup$shuffled.id[match(x, shuffle.gene.lookup$gene.id)])

  replacement.genes <- as.character(replacement.genes) # otherwise returns a factor or strings,
  # depending on versions
  df$locus_name[df$accession=='Ro18'] <- replacement.genes

  return(df)
}

plot_registration_for_exemplar_genes <- function(all.rep.shifted.data, GoIs) {
  # cut down to genes to plot
  tmp <- all.rep.shifted.data[all.rep.shifted.data$locus_name %in% GoIs]

  p.seperate <- ggplot(tmp, aes(x=timepoint, y=mean.cpm, color=locus_name, fill=locus_name,
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
  p.shifted <- ggplot(tmp, aes(x=shifted_time, y=mean.cpm,
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

  p.grid <- plot_grid(p.seperate, p.shifted, ncol=1, rel_heights = c(0.6, 0.4))
  # p.grid
  # ggsave('./test.pdf',
  #        plot=p.grid,
  #        width=3.2, height=6)

  return(p.grid)
}


make_data_heatmaps <- function(D.mean, D.scaled, D.registered, D.scaled.NR, D.scaled.R, D.registered.R) {

  p.mean <- make_heatmap(D.mean, 'mean expression')
  p.scaled <- make_heatmap(D.scaled, 'scaled mean expression')
  p.registered <- make_heatmap(D.registered, 'registered & scaled mean expression')

  p.scaled.NR <- make_heatmap(D.scaled.onlyNR, 'scaled mean expression (only not registered genes)')
  p.scaled.R <- make_heatmap(D.scaled.onlyR, 'scaled mean expression (only registered genes)')
  p.registered.R <- make_heatmap(D.registered.R, 'registered & scaled mean expression (only  registered genes)')


  p.all <- plot_grid(p.mean, p.scaled.NR, p.scaled, p.scaled.R, p.registered, p.registered.R, ncol=2)

  return(p.all)
}



make_heatmap_all <- function(D, title) {
  D$x.sample <- factor(D$x.sample, levels=unique(sort(D$x.sample)))
  D$y.sample <- factor(D$y.sample, levels=unique(sort(D$y.sample)))

  p <- ggplot(D, aes(x=x.sample, y=y.sample, fill=log(distance)))+
    geom_tile()+
    #scale_fill_viridis()+
    theme_classic()+
    facet_wrap(~title, ncol=1, scales='free')+
    theme(axis.text.x = element_text(angle=90,
                                     size=6),
          axis.text.y = element_text(size=6),
          plot.title = element_text(hjust=0.5, size=10),
          legend.position = 'top',
          legend.justification="right",
          legend.margin=margin(0,0,0,0),
          legend.box.margin = margin(0,0,-10,-10),
          legend.text=element_text(size=4, vjust=-0.5),
          legend.title = element_text(size=8),
          legend.key.height = unit(0.2, 'cm'),
    )+
    guides(fill=guide_colorbar(label.position='top'))+
    scale_fill_viridis()+
    ylab('')+
    xlab('')+
    ggtitle(title)

  return(p)
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

calculate_pairwise_sample_distance_simple <- function(dt) {
  # dt is the wide format expression of the two samples to be compared.
  # only genes which have data in both samples are considered (only
  # relevant for registration case).

  # filter to genes with data in both.
  dt <- na.omit(dt)
  # calculate distance
  dt$sq.diff <- (dt[, 1] - dt[, 2])^2

  d <- mean(dt$sq.diff)
  return(d)
}


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
  out.df <- data.table(data.frame('x.sample'=i.cols, 'y.sample'=j.cols, 'distance'=ds))
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

#brassica_name <- 'ro18_chiifu_apex'
get_data_all_symbol <- function(brassica_name) {

  rds_k <- 'klepikova' # 'klepikova

  # ------- changed by Ruth
  rdsdir <- 'final_data/rds/'

  rdspath_k <- paste0(rdsdir, rds_k, '.rds')
  rdspath_b <- paste0(rdsdir, brassica_name, '.rds')

  expression_b <- readRDS(rdspath_b)
  expression_k <- readRDS(rdspath_k)

  # id data
  # ------- changed by Ruth
  ID_TABLE <- readRDS(paste0('reference_data/ID_TABLE_brapa-v3.rds'))

  ID_TABLE <- unique(ID_TABLE[, c('locus_name', 'symbol', 'CDS.model')])

  # this screws up the number of copies of each gene unless reorder based on symbol name first!
  #ID_TABLE <- ID_TABLE[order(ID_TABLE$symbol, decreasing=T), ]
  #ID_TABLE <- ID_TABLE[!duplicated(ID_TABLE[c('locus_name', 'CDS.model')]), ]

  # add ATG locus info and cut down
  little.IDT <- unique(ID_TABLE[, c('CDS.model', 'locus_name')]) # has to be unique, because ID_TABLE has duplicate rows for multiple symbols!!
  expression_b <- merge(expression_b, little.IDT, by='CDS.model')

  #little.IDT <- ID_TABLE[, c('locus_name', 'symbol')]
  #names(little.IDT) <- c('CDS.model', 'symbol')
  #expression_k <- unique(merge(expression_k, little.IDT, by='CDS.model', allow.cartesian = T))
  expression_k$locus_name <- expression_k$CDS.model

  # cut down to only have genes with ATG locus present in both datasets
  common_symbols <- intersect(expression_k$locus_name, expression_b$locus_name)
  expression_k <- expression_k[expression_k$locus_name %in% common_symbols, ]
  expression_b <- expression_b[expression_b$locus_name %in% common_symbols, ]

  # join the two datasets into 1 & housekeeping
  expression <- rbind(expression_b, expression_k)

  # cut down to remove the 'blank' symbol
  expression <- expression[expression$locus_name!='',]

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

shorten_groups <- function(exp) {
  # get reps for klepikova and for brassica data
  exp <- data.table(exp)
  B <- exp[exp$accession != 'Col0']
  B[, c('j1', 'j2', 'j3', 'j4', 'j5', 'j6', 'rep'):=tstrsplit(sample_id, split='_')]
  B$rep[is.na(B$rep)] <- 1
  B[, c('j1', 'j2', 'j3', 'j4', 'j5', 'j6')] <- NULL

  A <- exp[exp$accession == 'Col0']
  A[, c('j1', 'j2','rep'):=tstrsplit(dataset, split='_')]
  A[, c('j1', 'j2')] <- NULL

  exp <- rbind(B,A)

  exp$ds <- mapvalues(exp$rep, from=c('1', '2', '3', '4'),
                      to=c('a', 'b', 'c', 'd'))
  exp <- data.table(exp)
  exp[, group:=paste(accession, sprintf('%02d', timepoint), ds, sep='-')]
  return(exp)
}
#models_o_I <- 'filt_models'

get_expression_oI <- function(rds_file, curr_GoIs, sumBrassicas) {

  # load rds and arabidopsis gene expression data into single df.
  master_exp <- get_data_all_symbol(rds_file)
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
    exp <- aggregate(norm.cpm~sample_id+accession+tissue+timepoint+dataset+group+locus_name, data=exp, sum)
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

# sample.id.cols <- c('accession','delta_time')
# gene.col <- c('locus_name')
# expression.col <- 'mean.cpm'
# dt <- mean_df
reformat_for_distance_calculation <- function(dt, sample.id.cols, gene.col, expression.col) {

  # concatenate sample.id columns to generate sample ids
  dt$sample.id <- dt[[sample.id.cols[1]]]
  if (length(sample.id.cols) > 1) {
    for (i in 2:length(sample.id.cols)) {
      # pad timepoint to 2 figures to help ordering
      if (class(dt[[sample.id.cols[i]]]) %in%  c('integer', 'numeric')) {
        dt[[sample.id.cols[i]]] <- str_pad(dt[[sample.id.cols[i]]], 2, pad='0')
      }

      dt$sample.id <- paste0(dt[['sample.id']], '-', dt[[sample.id.cols[i]]])
    }
  }

  # subset to just the relevant columns
  dt <- subset(dt, select=c('sample.id', gene.col, expression.col))

  # convert to wide format
  dt.w <- dcast(dt, locus_name ~ sample.id, value.var=eval(expression.col))

  return(dt.w)
}


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

  D <- dcast(D, group~locus_name, value.var='expression')

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
  #M <- na.omit(as.matrix(D[, -1]))
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
  apply(sc.M, 2, sd)

  #perp <- 1

  # ensure graph file exists
  if (!dir.exists(save.dir)) {
    dir.create(save.dir)
  }

  perp <- 1
  for (perp in c(1, 2, 3, 4, 5, 6, 8)) {
    # set to 8000 iterations -≥ pretty consistently minimises the error to ~0.35
    tsne <- Rtsne(sc.M, dims=2, perplexity=perp, verbose=T, max_iter=16000, pca=F)
    p1 <- tSNEPlot(tsne, Labels) + guides(color=F)
    tsne <- Rtsne(sc.M, dims=2, perplexity=perp, verbose=T, max_iter=16000, pca=F)
    p2 <- tSNEPlot(tsne, Labels) + guides(color=F)
    tsne <- Rtsne(sc.M, dims=2, perplexity=perp, verbose=T, max_iter=16000, pca=F)
    p3 <- tSNEPlot(tsne, Labels) + guides(color=F)
    tsne <- Rtsne(sc.M, dims=2, perplexity=perp, verbose=T, max_iter=16000, pca=F)
    p4 <- tSNEPlot(tsne, Labels) + guides(color=F)
    p <- ggarrange(p1, p2, p3, p4,
                   labels=c('a', 'b', 'c', 'd'),
                   ncol=2, nrow=2)
    p
    ggsave(paste0(save.dir, '/', plot_name, '_p=', perp,'.pdf'), scale=1.5)
  }

  # try with PCA before tSNE, as it seems to give an ok gradient timecourse.
  # for (perp in c(2,3,4,5,68,12,15)) {
  #   # set to 8000 iterations -≥ pretty consistently minimises the error to ~0.35
  #   tsne <- Rtsne(sc.M, dims=2, perplexity=perp, verbose=T, max_iter=16000, pca=T, initial_dims=12)
  #   p1 <- tSNEPlot(tsne, Labels) + guides(color=F)
  #   tsne <- Rtsne(sc.M, dims=2, perplexity=perp, verbose=T, max_iter=16000, pca=T, initial_dims=12)
  #   p2 <- tSNEPlot(tsne, Labels) + guides(color=F)
  #   tsne <- Rtsne(sc.M, dims=2, perplexity=perp, verbose=T, max_iter=16000, pca=T, initial_dims=12)
  #   p3 <- tSNEPlot(tsne, Labels) + guides(color=F)
  #   tsne <- Rtsne(sc.M, dims=2, perplexity=perp, verbose=T, max_iter=16000, pca=T, initial_dims=12)
  #   p4 <- tSNEPlot(tsne, Labels) + guides(color=F)
  #   p <- ggarrange(p1, p2, p3, p4,
  #                  labels=c('a', 'b', 'c', 'd'),
  #                  ncol=2, nrow=2)
  #   ggsave(paste0('../graphs/', rds_file, '/tSNE_plots/', plot_name, '_p=', perp,'_PCA.pdf'), scale=1.5)
  # }
}

tSNEPlot <- function(tsne, Labels) {
  position <- data.frame(tsne$Y)
  position <- cbind(position, Labels)
  names(position) <- c('x', 'y', 'label')

  position <- data.table(position)
  position[, c('accession', 'timepoint', 'letter'):=tstrsplit(label, split='-')]
  p <- ggplot(position, aes(x=x, y=y, color=accession, label=paste(timepoint, letter, sep='-')))+
    geom_point()+
    geom_text_repel(size=2)
  return(p)
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



impute_arabidopsis_values <- function(shifted.mean_df) {
  # Arabidopsis gene expression profiles are shifted all over. Need to impute times at set of common timepoints
  # in order to allow sample distance comparison to Ro18.
  #
  # Ro18 genes haven't been shifted around, therefore imputed timepoints are realtive to Ro18 timepoints.
  # We ONLY have ro18 observations for 11, 13, 15, ... 35
  #
  # Col0 can have had different shifts applied to them, so need to impute them to a common scale, and to compare to Ro18
  #
  # therefore, we only need to impute Arabidopsis gene expression, which will be compared to the Ro18 timepoints.
  # BUT don't want to discard any col0 data, so want to generate imputed time observations for col0 from minimum
  # to maximum shifted timepoints for Col0 not just for 11,13,15, etc RO18 observations.


  # sanity plotting - a registered one
  # ggplot(shifted.mean_df[shifted.mean_df$locus_name=='BRAA01G001090.3C', ], aes(x=shifted_time, y=mean.cpm, color=accession))+
  #   geom_point()
  # # an unregistered one
  # ggplot(shifted.mean_df[shifted.mean_df$locus_name=='BRAA01G003140.3C', ], aes(x=shifted_time, y=mean.cpm, color=accession))+
  #   geom_point()

  # The imputed col0 times going to extimate gene expression for
  imputed.timepoints <- round(seq(min(shifted.mean_df$shifted_time), max(shifted.mean_df$shifted_time)))

  out.list <- list()
  out.list <- c(out.list, list(shifted.mean_df[shifted.mean_df$accession=='Ro18']))

  curr.gene <- 'BRAA01G001090.3C' #unique(shifted.mean_df$locus_name)[17]
  count <- 0
  for (curr.gene in unique(shifted.mean_df$locus_name)) {
    if (count %% 100 == 0) {
      print(paste0(count, ' / ', length(unique(shifted.mean_df$locus_name))))
    }

    # get the current gene expression data
    curr.df <- shifted.mean_df[shifted.mean_df$locus_name==curr.gene, ]

    # skip over this one, if not assigned a best shift value, because brassica gene not expressed
    # if (all(is.na(curr.df$mean.cpm))) {
    #   next
    # }

    ara.df <- curr.df[curr.df$accession=='Col0',]
    bra.df <- curr.df[curr.df$accession=='Ro18',]

    interp.ara.df <- data.table(data.frame('locus_name'=curr.gene, 'accession'='Col0', 'tissue'='apex', 'timepoint'=NA,
                                           'stretched.time.delta'= NA, 'shifted_time'=imputed.timepoints,
                                           'is.registered'= unique(ara.df$is.registered)[1]))

    # for each brassica timepoint, interpolate the comparible arabidopsis expression
    # by linear interpolation between the neighbouring 2 ara values. If not between 2 ara values
    # because shifted outside comparible range, set to NA
    interp.ara.df$mean.cpm <- sapply(imputed.timepoints, interpolate_brassica_comparison_expression, bra.dt=ara.df)


    # sanity testing - line is interpolated.
    # ggplot(curr.df, aes(x=shifted_time, y=mean.cpm, color=accession))+
    #   geom_point()+
    #   geom_line(data=interp.ara.df)

    out.list <- c(out.list, list(interp.ara.df))
    count <- count+1
  }
  out.df <- do.call('rbind', out.list)
  return(out.df)
}

# test <- mean_df
# stretch_factor
# min_num_overlapping_points
# shift_extreme
get_extreme_shifts_for_all <- function(test, stretch_factor, min_num_overlapping_points, shift_extreme) {
  # wrapper for calc_extreme_shifts to be able to move it out of the loop so don't calculate for every gene.

  #min_num_overlapping_points <- 5 # bound the extreme allowed shifts, such than at least this many timepoints are being compared for both accessions.
  # cut data.table to a single gene
  curr_sym <- unique(test$locus_name)[1]
  test <- test[test$locus_name==curr_sym, ]

  # transform timepoint to be time from first timepoint
  test[, delta_time:=timepoint - min(timepoint), by=.(accession)]
  # apply stretch_factor to the arabidopsis, leave the rapa as is
  test$delta_time[test$accession=='Col0'] <- test$delta_time[test$accession=='Col0']*stretch_factor

  # calculate min shift and max time shift, which still allows overlap of at least 5 times to be compared from whichever accession will be considering fewer timepoints from.
  # Shift is applied to the arabidopsis - so the 5th largest arabidopsis time is the biggest -ve shift can be applied
  # and the biggest shift which can be applied is to make the 5th smallest arabidopsis time == largest brassica time
  #setorder(test, delta_time)
  #min_shift <- min(test$delta_time[test$accession=='Ro18']) - test$delta_time[test$accession=='Col0'][length(test$delta_time[test$accession=='Col0'])-4]
  #max_shift <- max(test$delta_time[test$accession=='Ro18']) - test$delta_time[test$accession=='Col0'][5]

  M <- calc_extreme_shifts(test, min_num_overlapping_points, shift_extreme)
  return(M)
}



# for all genes, in mean_df, get the scores, and the shifts which result in them after
# stretching arabidopsis by "stretch_factor, and applying shifts forward and
# backward, with the extremes defined by allowing 5 overlapping points for comparison.
# mean_df <- to.shift.df
# stretch_factor <- 2
# do_rescale=F
calculate_all_best_shifts <- function(mean_df, stretch_factor, do_rescale, min_num_overlapping_points, shift_extreme) {
  symbols <- c()
  num_points <- c()
  #curr_sym <- 'TT16'
  all_scores_list <- rep(list(0), length(unique(mean_df$locus_name)))
  length(unique(mean_df$locus_name))


  # get the extreme shifts which can be applied to the genes
  M <- get_extreme_shifts_for_all(mean_df, stretch_factor, min_num_overlapping_points, shift_extreme)
  min_shift <- M[[1]]
  max_shift <- M[[2]]

  count <- 0
  curr_sym <- unique(mean_df$locus_name)[2]

  i = 1
  for (i in 1:length(unique(mean_df$locus_name))) {
  #for (curr_sym in unique(mean_df$locus_name)) {
    curr_sym <- unique(mean_df$locus_name)[i]
    if (count %% 100 == 0) {
      print(paste0(count, ' / ', length(unique(mean_df$locus_name))))
    }

    # out is mean SSD between arabidopsis, and interpolated brassica (interpolated between 2 nearest points)
    # ggplot(mean_df[mean_df$locus_name==curr_sym,], aes(x=timepoint, y=mean.cpm, color=accession))+
    #   geom_point()

    ### get "score" for all the candidate shifts - score is mean error / brassica expression for compared points.
    ### if timepoints don't line up, brassica value is linearly imputed
    out <- get_best_shift_new(curr_sym, mean_df, stretch_factor, do_rescale, min_shift, max_shift, testing=FALSE)

    best_shift <- out$shift[out$score==min(out$score)]
    if (length(best_shift) > 1) {
      if (max(out$score)=='Inf') { # can get inf score if brassica gene note expressed in the comparison
        next
      } else {
        # if ties for the best shift applied, apply the smaller absolute one
        best_shift <- best_shift[abs(best_shift) == min(abs(best_shift))]
      }
    }

    all_scores <- out
    all_scores_list[[i]] <- all_scores
    symbols <- c(symbols, curr_sym)
    #num_points <- c(num_points, best_shift)
    #print(best_shift)

    count <- count + 1
  }
  #shift_results <- data.frame(symbol=symbols, num.points = num_points)
  all_scores_df <- do.call('rbind', all_scores_list)

  return(all_scores_df)
}

#test <- mean_df
get_best_shift <- function(curr_sym, test) {
  # here statistic used to assess best is 1 / mean(squared difference)

  min_shift <- 5

  # cut to get a single symbol
  test <- test[test$locus_name==curr_sym, ]

  # get the mean expression vectors for the current gene
  AtVec <- test$mean.cpm[test$accession=='Col0']
  BrVec <- test$mean.cpm[test$accession=='Ro18']

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

# wroking on new implementation of the funciton
# curr_sym <- 'BRAA02G015410.3C'
#test <- mean_df
# stretch_factor <- 1
# do_rescale <- TRUE
#testing=T

get_best_shift_new <- function(curr_sym, test, stretch_factor, do_rescale, min_shift, max_shift, testing=FALSE) {
  # for the current gene, and current stretch_factor, calculate the score for all
  # shifts, and return the scores for all as a table, and the value of the optimal shift.

  # Shift extremes are defined s.t. at least 5 points are compared.

  # do_rescale == TRUE, means apply "scale" to compared points for each shift. ==FALSE, means use original mean expression data

  num.shifts <- 10 # the number of different shifts to be considered.


  test <- test[test$locus_name==curr_sym, ]

  # transform timepoint to be time from first timepoint
  test[, delta_time:=timepoint - min(timepoint), by=.(accession)]
  # apply stretch_factor to the arabidopsis, leave the rapa as is
  test$delta_time[test$accession=='Col0'] <- test$delta_time[test$accession=='Col0']*stretch_factor

  all_scores <- rep(0, num.shifts)
  all.ara.mean <- rep(0, num.shifts)
  all.bra.mean <- rep(0, num.shifts)
  all.ara.sd <- rep(0, num.shifts)
  all.bra.sd <- rep(0, num.shifts)

  all.shifts <- seq(min_shift, max_shift, length.out=num.shifts)
  if (!(0 %in% all.shifts)) {
    all.shifts <- c(all.shifts, 0) # include 0 shift in candidates.
  }
  i=1
  for (i in 1:length(all.shifts)) {

    curr.shift <- all.shifts[i]

    #print('line 1676')
    #print(curr.shift)

    # shift the arabidopsis expression timeings
    test$shifted_time <- test$delta_time
    test$shifted_time[test$accession=='Col0'] <- test$delta_time[test$accession=='Col0'] + curr.shift

    #### test plot - of shifted, UNNORMALISED gene expression
    # if (testing==TRUE) {
    #   p <- ggplot(test, aes(x=shifted_time, y=mean.cpm, color=accession))+
    #     geom_point()+
    #     ggtitle(paste0('shift : ', curr.shift))
    #   p
    #   ggsave(paste0('./testing/', curr.shift, '.pdf'))
    # }


    # cut down to just the arabidopsis and brassica timepoints which compared
    test <- get_compared_timepoints(test)
    compared <- test[test$is.compared==TRUE, ]

    # renormalise expression using just these timepoints?
    if (do_rescale==TRUE) {
      # record the mean and sd of the compared points, used for rescaling
      # in "apply shift" function
      ara.mean <- mean(compared$mean.cpm[compared$accession=='Col0'])
      bra.mean <- mean(compared$mean.cpm[compared$accession=='Ro18'])
      ara.sd <- sd(compared$mean.cpm[compared$accession=='Col0'])
      bra.sd<- sd(compared$mean.cpm[compared$accession=='Ro18'])

      # do the transformation for here
      if ((ara.sd != 0) & (bra.sd != 0)) { # if neither are 0, so won't be dividing by 0 (which gives NaNs)
        compared[, mean.cpm:=scale(mean.cpm, scale=TRUE, center=TRUE), by=.(accession)]
      } else { # if at least one of them is all 0
        ara.compared <- compared[compared$accession=='Col0',]
        bra.compared <- compared[compared$accession=='Ro18',]
        if((ara.sd == 0) & (bra.sd != 0)) { # if only ara.sd==0
          bra.compared[, mean.cpm:=scale(mean.cpm, scale=TRUE, center=TRUE), by=.(accession)]
        }
        if ((ara.sd != 0) & (bra.sd == 0)) { # if only bra.sd == 0
          ara.compared[, mean.cpm:=scale(mean.cpm, scale=TRUE, center=TRUE), by=.(accession)]
        }
        # if both are all 0, then do nothing.
        compared <- rbind(ara.compared, bra.compared)
      }

    } else {
      # if didn't rescale expression for comparison, record values s.t. (x - xmean) / x_sd = x
      ara.mean <- 0
      bra.mean <- 0
      ara.sd <- 1
      bra.sd<- 1
    }

    ### test plot of shifted, and normalised gene expression
    if (testing==TRUE) {
      p <- ggplot(compared, aes(x=shifted_time, y=mean.cpm, color=accession))+
        geom_point()+
        ggtitle(paste0('shift : ', curr.shift))
        ggsave(paste0('./testing/',stretch_factor, '-', curr.shift, '.pdf'))
    }

    # for each arabidopsis timepoint, linear interpolate between the two nearest brassica timepoints
    ara.compared <- compared[compared$accession=='Col0']
    bra.compared <- compared[compared$accession=='Ro18']

    ara.compared$pred.bra.expression <- sapply(ara.compared$shifted_time, interpolate_brassica_comparison_expression, bra.dt=bra.compared)

    # calculate the score, using the (interpolated) predicted.bra.expression, and the observed arabidopsis expression
    # score = mean ((observed - expected)**2 )
    score <- calc_score(ara.compared$mean.cpm, ara.compared$pred.bra.expression)

    if (is.na(score)) {
      print('error in get_best_shift_new(): got a score of NA for gene:')
      print(curr_sym)
      print(paste0('with curr.shift=', curr.shift))
      stop()
    }

    all_scores[i] <- score
    all.ara.mean[i] <- ara.mean
    all.bra.mean[i] <- bra.mean
    all.ara.sd[i] <- ara.sd
    all.bra.sd[i] <- bra.sd
  }

  out <- data.table(data.frame('gene'=curr_sym, 'stretch'=stretch_factor, 'shift'=all.shifts, 'score'=all_scores,
                               'ara.compared.mean'=all.ara.mean, 'bra.compared.mean'=all.bra.mean,
                               'ara.compared.sd'=all.ara.sd, 'bra.compared.sd'=all.bra.sd))
  return(out)
}


calc_num_overlapping_points <- function(shift, original) {
  # calculate the number of overlapping points for the species with the fewer overlapping points if the current "shift" is
  # applied to the col0 delta timepoints.
  original$shifted_time[original$accession=='Col0'] <- original$delta_time[original$accession=='Col0'] + shift
  original <- get_compared_timepoints(original)
  original[, num.compared:=sum(is.compared), by=.(accession)]

  return(min(original$num.compared))
}

calc_extreme_shifts <- function(test, min_num_overlapping_points, shift_extreme) {
  # calculate the minimum and maximum shifts can apply to Col-0 after the stretch transformation, whilst
  # preserving the criteria that at least min_num_overlapping_points are being compared from both accessions.

  original <- copy(test)
  original$shifted_time <- original$delta_time

  # print('line 1803')
  # print(original)

  # -ve extreme shift will be -1*exactly the difference between 1 of the stretched Col0 timepoints, and the smallest ro18 timepoint
  # +ve extreme will be the difference between 1 of the col0 timepoints, and the maximum Ro18 timepoint
  neg_extreme_candidate <- -1*(original$delta_time[original$accession=='Col0'] - min(original$delta_time[original$accession=='Ro18']))
  pos_extreme_candidates <- max(original$delta_time[original$accession=='Ro18']) - original$delta_time[original$accession=='Col0']

  # of these candidates, find the most extreme values which mainting the required number of overlapping timepoints to be considered.
  num_overlapping_points <- sapply(neg_extreme_candidate, FUN=calc_num_overlapping_points, original=original)
  if (all(num_overlapping_points < min_num_overlapping_points)) {
    stop(paste0('calc_extreme_shifts():\nafter applying stretch factor:', stretch, ' to ', transformed.timecourse, ', none of the considered shifts have ',
                 'min_num_overlapping_points (', min_num_overlapping_points, ') overlapping timepoints with the other timecourse!\n',
                "maybe try a smaller stretch, and double check you're applying it to the correct timecourse." ))
  }

  neg.extreme <- min(neg_extreme_candidate[num_overlapping_points >= min_num_overlapping_points])

  num_overlapping_points <- sapply(pos_extreme_candidates, FUN=calc_num_overlapping_points, original=original)
  pos_extreme <- max(pos_extreme_candidates[num_overlapping_points >= min_num_overlapping_points])

  # hard code maximum and minimum allowed shifts, as noticed spurious registrations when too extreme shifts
  # allowed
  if (neg.extreme < (-1*shift_extreme)) {
    neg.extreme <- -1 * shift_extreme
  }
  if (pos_extreme > 1*shift_extreme) {
    pos_extreme <- shift_extreme
  }

  return(list(neg.extreme, pos_extreme))
}



apply_stretch <- function(mean_df, best_shifts) {
  # gets the applied stretch from the best_shifts df
  test <- copy(mean_df)

  # get the stretch factor applied in the best_shifts
  #stopifnot(length(unique(best_shifts$stretch))==1)
  #stretch_factor <- unique(best_shifts$stretch)

  #stopifnot(length(unique(mean_df$locus_name))==length(best_shifts$gene)) # best_shifts should have 1 row per gene
  # if (length(unique(mean_df$locus_name)) != length(best_shifts$gene)) {
  #   print('ERROR!:')
  #   print(length(unique(mean_df$locus_name)))
  #   print(length(best_shifts$gene))
  #   stop()
  # }


  # stretch the arabidopsis expression data, leave the rapa as is
  test[, delta_time:=timepoint - min(timepoint), by=.(accession)]
  ro18.test <- test[test$accession=='Ro18',]
  col0.test <- test[test$accession=='Col0']
  #test$delta_time[test$accession=='Col0'] <- test$delta_time[test$accession=='Col0']*stretch_factor
  col0.test <- merge(col0.test, best_shifts[, c('gene', 'stretch')], by.x='locus_name', by.y='gene')
  col0.test$delta_time <- col0.test$delta_time * col0.test$stretch
  col0.test$stretch <- NULL
  test <- rbind(ro18.test, col0.test)

  # record the stretched times (before indiv shifting applied)
  test$stretched.time.delta <- test$delta_time # record the time (from start of timecourse) after stretching,
  test$shifted_time <- test$delta_time
  # after stretching, add the time to the first datapoint (7d for ara, 11d for ro18) back on
  test$shifted_time[test$accession=='Col0'] <- test$shifted_time[test$accession=='Col0'] + 11 #7
  test$shifted_time[test$accession=='Ro18'] <- test$shifted_time[test$accession=='Ro18'] + 11
  test$delta_time <- NULL

  return(test)
}


# mean_df <- all.data.df
# best_shifts
apply_best_shift <- function(mean_df, best_shifts) {
  # take unregistered expression over time, and the best shifts, and
  # return the registered expression over time for each gene

  test <- copy(mean_df)

  test <- apply_stretch(mean_df, best_shifts)

  # normalise the expression data (if was normalised when calculating the expression data, is recorder in the
  # .compared.mean, and .compared.sd columns. If no normalisation was carried out, then these should have values of 0,
  # and 1).
  if (!(all(unique(best_shifts$ara.compared.mean) == 0)) |
      !(all(unique(best_shifts$bra.compared.mean) == 0))) {
    print('Normalising expression by mean and sd of compared values...')
    test <- apply_best_normalisation(test, best_shifts)
    print('done!')

  } else { # if no scaling carried out DURING the registration step
    print('No normalisation was carried out DURING registration (though may have been, prior to Col-Ro18 comparison)')
    test <- test
  }

  print('applying best shift...')

  # for each gene, shift the arabidopsis expression by the optimal shift found previously
  curr.gene <- 'BRAA01G000040.3C'
  #curr.gene <- unique(test$locus_name)[1]
  for (curr.gene in unique(test$locus_name)) {
    #print(curr.gene)
    curr.best.shift <- best_shifts$shift[best_shifts$gene==curr.gene]
    test$shifted_time[test$accession=='Col0' & test$locus_name==curr.gene] <- test$shifted_time[test$accession=='Col0' & test$locus_name==curr.gene] + curr.best.shift

    # tmp <- test[test$locus_name==curr.gene]
    # ggplot(tmp, aes(x=shifted_time, y=mean.cpm, color=accession))+
    #    geom_point()
  }



  print('done!')

  return(test)
}


apply_best_normalisation <- function(test, best_shifts) {
  # for each gene, in each accession (Ro18 and COl0) normalise by the mean and standard deviation of the compared points.
  # if the gene wasn't compared, set the expresion value to NA

  # curr.gene <- 'BRAA01G001470.3C'
  count = 0
  for (curr.gene in unique(test$locus_name)) {
    if (count %% 100 == 0) {
      print(paste0(count, ' / ', length(unique(test$locus_name))))
    }
    ara.mean <- best_shifts$ara.compared.mean[best_shifts$gene==curr.gene]
    bra.mean <- best_shifts$bra.compared.mean[best_shifts$gene==curr.gene]
    ara.sd<- best_shifts$ara.compared.sd[best_shifts$gene==curr.gene]
    bra.sd<- best_shifts$bra.compared.sd[best_shifts$gene==curr.gene]

    # if was compared
    if (length(ara.mean) != 0) {
      if (ara.sd != 0) { # don't want to divide by 0
        test$mean.cpm[test$locus_name==curr.gene & test$accession=='Col0'] <- (test$mean.cpm[test$locus_name==curr.gene & test$accession=='Col0'] - ara.mean) / ara.sd
      } else {
        test$mean.cpm[test$locus_name==curr.gene & test$accession=='Col0'] <- (test$mean.cpm[test$locus_name==curr.gene & test$accession=='Col0'] - ara.mean)
      }

      if (bra.sd !=0) { # don't want to divide by 0
        test$mean.cpm[test$locus_name==curr.gene & test$accession=='Ro18'] <- (test$mean.cpm[test$locus_name==curr.gene & test$accession=='Ro18'] - bra.mean) / bra.sd
      } else {
        test$mean.cpm[test$locus_name==curr.gene & test$accession=='Ro18'] <- (test$mean.cpm[test$locus_name==curr.gene & test$accession=='Ro18'] - bra.mean)
      }
      if (any(is.na(test$mean.cpm))) {
        print('have NAs in mean.cpm after rescaling in apply best_normalisation() for gene :')
        print(unique(test$locus_name))
        stop()
      }
    } else {
      test$mean.cpm[test$locus_name==curr.gene & test$accession=='Col0'] <- NA
      test$mean.cpm[test$locus_name==curr.gene & test$accession=='Ro18'] <- NA
    }
    count <- count + 1
  }

  return(test)
}

# my.scale <- function(V) {
#   V <- (V - mean(V)) / sd(V)
#   return(V)
# }

# test <- curr.data.df
get_compared_timepoints <- function(test) {
  # flag the arabidopsis timepoints which overlap the brassica timecourse, and so will be compared
  bra.min <- min(test$shifted_time[test$accession=='Ro18'])
  bra.max <- max(test$shifted_time[test$accession=='Ro18'])

  # get the arabidopsis times which used
  test$is.compared <- FALSE
  test$is.compared[(test$accession=='Col0' & (test$shifted_time >= bra.min & test$shifted_time <=bra.max))] <- TRUE

  # get the extreme brassica times which used - bigger or equal than Ara max, and smaller or equal than Ara min, because have to project
  #  Ara onto Bra
  ara.max <- max(test$shifted_time[test$accession=='Col0' & test$is.compared==TRUE])
  ara.min <- min(test$shifted_time[test$accession=='Col0' & test$is.compared==TRUE])
  bra.max <- max_is_compared_to_arabidopsis(ara.max, test[test$accession=='Ro18', ])
  bra.min <- min_is_compared_to_arabidopsis(ara.min, test[test$accession=='Ro18', ])

  # use these to get all the brassica times which used
  test$is.compared[(test$accession=='Ro18' & (test$shifted_time >= bra.min & test$shifted_time <=bra.max))] <- TRUE

  return(test)
}

#ara.expression <- ara.compared$mean.cpm
#bra.expression <- ara.compared$pred.bra.expression
calc_score <- function(ara.expression, bra.expression) {
  # (sum(observed-expected)**2)
  # take mean, because going to be comparing variable number of datapoints.
  # if don't regularise / penalise for shift applied, then like uniform prior on it.
  # maybe should penalise for comparing fewer timepoints?
  # divide by the ara.expression as filtered already to make sure expressed

  d <- (ara.expression - bra.expression)**2 #/ abs(ara.expression)
  return(mean(d))
}

# arabidopsis.time <- ara.max
#arabidopsis.time <- ara.max
#bra.dt <- test[test$accession=='Ro18',]
max_is_compared_to_arabidopsis <- function(arabidopsis.time, bra.dt) {
  # return the largest brassica time which is used in comparison t the arabidopsis time
  # the smallest one greater to or equal to arabidopsis time

  # if using for rep data, then repeats of the same points screws it up
  bra.dt <- unique(subset(bra.dt, select=c('timepoint', 'shifted_time')))

  # return the bra dt shifted timepoint which is greater than, or equal to the Ara time.
  bra.dt$diff <- bra.dt$shifted_time - arabidopsis.time
  candidates <- bra.dt[bra.dt$diff>=0,]
  bra.max.time <- candidates$shifted_time[candidates$diff==min(candidates$diff)]

  return(bra.max.time)

  #setorder(bra.dt, diff)
  #nearest.points <- bra.dt[1:2,]
  #setorder(nearest.points, shifted_time)
  #return(nearest.points$shifted_time[2])
}

#arabidopsis.time=ara.min
min_is_compared_to_arabidopsis <- function(arabidopsis.time, bra.dt) {
  # return the smallest brassica time which is used in comparison t the arabidopsis time
  # the biggest one smaller than or equal to the arabidopsis time

  # if using for rep data, then repeats of the same points screws it up
  bra.dt <- unique(subset(bra.dt, select=c('timepoint', 'shifted_time')))

  bra.dt$diff <- bra.dt$shifted_time - arabidopsis.time
  candidates <- bra.dt[bra.dt$diff<=0,]
  bra.min.time <- candidates$shifted_time[candidates$diff==max(candidates$diff)]
  return(bra.min.time)

}

# arabidopsis.time <- 2
# bra.dt <- ara.df
interpolate_brassica_comparison_expression <- function(arabidopsis.time, bra.dt) {

  # arabidopsis time is outside of the range of the bra.dt shifted timepoints
  bra.dt$diff <- bra.dt$shifted_time - arabidopsis.time

  # if outside of comparible range (time is smaller than all bra.dt time or bigger than all)
  if (all(bra.dt$diff > 0) | all(bra.dt$diff < 0) ) {
    return(NA)
  }

  # otherwise,  cut down brassica observations to the two nearest timepoints to the arabidopsis time
  bra.dt$diff <- abs(bra.dt$shifted_time - arabidopsis.time)
  setorder(bra.dt, diff)
  nearest.points <- bra.dt[1:2,]

  # linearly interpolate between these points to estimate the comparison expression value
  setorder(nearest.points, shifted_time) # so [1] is earlier time
  time.diff <- nearest.points$shifted_time[2] - nearest.points$shifted_time[1] #
  expression.diff <- nearest.points$mean.cpm[2] - nearest.points$mean.cpm[1]
  grad <- expression.diff / time.diff
  pred.expression <- nearest.points$mean.cpm[1] + (nearest.points$diff[1]) * grad

  return(pred.expression)
}


#get_shifted_expression(shift_results, exp)
get_shifted_expression <- function(shift_results, exp) {
  cur_gene <- 'BRAA01G010430.3C'
  shifted_exp <- list()
  for (cur_gene in unique(shift_results$symbol)) {

    # cut to get a single symbol
    test <- exp[exp$locus_name==cur_gene, ]

    # ggplot(test, aes(x=timepoint, y=norm.cpm, color=accession))+
    #   geom_point()

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
    # ggplot(df, aes(x=shifted_time, y=norm.cpm, color=accession))+
    #   geom_point()

    #adf <- data.frame(accession='Col0', timepoint=at_times, sc.norm.cpm=AtVec, symbol=cur_gene, shift=num_points/3)
    #bdf <- data.frame(accession='Ro18', timepoint=br_times, sc.norm.cpm=BrVec, symbol=cur_gene, shift=num_points/3)
    #df <- rbind(adf, bdf)
    shifted_exp <- c(shifted_exp, list(df))
  }
  shifted_exp <- do.call('rbind', shifted_exp)

  return(shifted_exp)
}
