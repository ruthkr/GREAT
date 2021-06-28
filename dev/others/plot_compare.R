plot_goI_expression <- function(summed.GoIs.df) {

  # make plot of gene expression in Col0, and in Ro18.
  # truncate data so can see expression nicely on the same scale
  # summed.GoIs.df <- summed.GoIs.df[summed.GoIs.df$timepoint <=21,]

  summed.GoIs.df[, scaled.cpm:=my_scale(mean_cpm), by=.(Ara.id, accession)]
  morphology.equiv.df <- data.frame('accession'=c('Col-0', 'DH'), 'floral.transition.time'=c(14, 35))
  summed.GoIs.df$accession <- as.character(summed.GoIs.df$accession)
  summed.GoIs.df$accession[summed.GoIs.df$accession=='Col0'] <- 'Col-0'
  summed.GoIs.df$accession[summed.GoIs.df$accession=='Ro18'] <- 'DH'

  curr.acc <- 'Col-0'
  plot.list <- list()
  for (curr.acc in c('Col-0', 'DH')) {

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

plot_goI_expression_scaled <- function(summed.GoIs.df) {

  # make plot of gene expression in Col0, and in Ro18.
  # truncate data so can see expression nicely on the same scale
  # summed.GoIs.df <- summed.GoIs.df[summed.GoIs.df$timepoint <=21,]

  summed.GoIs.df[, scaled.cpm:=my_scale(mean_cpm), by=.(Ara.id, accession)]
  morphology.equiv.df <- data.frame('accession'=c('Col-0', 'DH'), 'floral.transition.time'=c(14, 35))
  summed.GoIs.df$accession <- as.character(summed.GoIs.df$accession)
  summed.GoIs.df$accession[summed.GoIs.df$accession=='Col0'] <- 'Col-0'
  summed.GoIs.df$accession[summed.GoIs.df$accession=='Ro18'] <- 'DH'

  curr.acc <- 'Col-0'
  plot.list <- list()
  for (curr.acc in c('Col-0', 'DH')) {

    curr.p <- ggplot(summed.GoIs.df[summed.GoIs.df$accession==curr.acc,],
                     aes(x=timepoint, y=scaled.cpm, color=Ara.name, fill=Ara.name))+
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
      curr.p <- curr.p + theme(axis.title.x=element_blank())
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
  AGL24.df <- registered.plot.df[registered.plot.df$locus_name=='AGL24'] #&
                                   # (registered.plot.df$shifted_time <=35) &
                                   # (registered.plot.df$shifted_time >=13),]
  AP1.df <- registered.plot.df[registered.plot.df$locus_name=='AP1'] #&
                                 # (registered.plot.df$shifted_time <=25) &
                                 # (registered.plot.df$shifted_time >=11),]
  AP3.df <- registered.plot.df[registered.plot.df$locus_name=='AP3'] #&
                                 # (registered.plot.df$shifted_time <=25) &
                                 # (registered.plot.df$shifted_time >=13),]
  LFY.df <- registered.plot.df[registered.plot.df$locus_name=='LFY'] #&
                                 # (registered.plot.df$shifted_time <=19) &
                                 # (registered.plot.df$shifted_time >=10),]
  SOC1.df <- registered.plot.df[registered.plot.df$locus_name=='SOC1'] #&
                                  # (registered.plot.df$shifted_time <=31) &
                                  # (registered.plot.df$shifted_time >=10),]
  SPL8.df <- registered.plot.df[registered.plot.df$locus_name=='SPL8']
  SVP.df <- registered.plot.df[registered.plot.df$locus_name=='SVP']
  SPL5.df <- registered.plot.df[registered.plot.df$locus_name=='SPL5']
  registered.plot.df <- rbind(AGL24.df, AP1.df, AP3.df, LFY.df, SOC1.df, SPL8.df, SVP.df, SPL5.df)

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
    # scale_color_brewer(palette = "Set2", direction = 1) +
    # scale_fill_brewer(palette = "Set2", direction = 1) +
    theme(legend.position = 'top',
          legend.title = element_blank(),
          axis.title = element_text(size=10),
          axis.text=element_text(size=6),
          strip.text=element_text(face='italic'),
          legend.margin=margin(21,0,0,0))

  return(p.registered)
}

plot_registered_GoIs_for_comparible_timepoints_new <- function(registered.plot.df) {

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
    # scale_color_brewer(palette = "Set2", direction = 1) +
    # scale_fill_brewer(palette = "Set2", direction = 1) +
    theme(legend.position = 'top',
          legend.title = element_blank(),
          axis.title = element_text(size=10),
          axis.text=element_text(size=6),
          strip.text=element_text(face='italic'),
          legend.margin=margin(21,0,0,0))

  return(p.registered)
}
