rm(list=ls())
#library(Rtsne)
# library(data.table)
#library(ggrepel)
# library(plyr)
#library(ggpubr)
# library(ggplot2)
# library(viridis)
# library(cowplot)
# library(stringr)
# library(splines)

## once have used ballgown_Arabidopsis_get_comparison_genes.R to make a .tsv files of candidate mapping genes,
# and the methods used to derive them, use this to make tSNE comparison plot.

# here, register the gene expression data for the real data, and for shuffled (randomised) data.

### IF RUNNING ON CLUSTER
# args <- commandArgs(trailingOnly = T)
# jobNum <- toString(args[1])
# do.initial.rescale <- toString(args[2])
# do.register.rescale <- toString(args[3])
# shuffle.type <- toString(args[4])
# outdir.string <- toString(args[5])
# setwd('/jic/research-projects/bravo/alex/rapa_paper/scripts/')
# source("./registration_functions.R")


### IF RUNNING LOCALLY
jobNum <- 1 # if running lots of times on cluster to get enough shuffled results, then jobNum is used to prevent output overwriting
setwd('/Volumes/Research-Projects/bravo/shared/to_Shannon/scripts/')
source("./chiifu_v3_scripts/registration_functions.R")
do.initial.rescale <- 'nope' # should be 'rescale' if want to use scaled df for registration, rather than mean_df
do.register.rescale <- 'rescale' # should be 'rescale' if want to rescale using only the overlapping points during
                              # registration
shuffle.type <- 'shuffle.expression' # whether should shuffle by shuffling the gene ids compared, or by shuffling the gene expression
                               # for each gene c('shuffle.genes', or 'shuffle.expression)
outdir.string <- 'TESTING_rescale_as_register___shuffled_g_v4__'





# calculate the scores associated with each candidate shift
# Arabidopsis observations = 7d -> 16d = 0d -> 9d
# Ro18 observations = 11d -> 35d = 0d -> 24d
# therefore if same start and end points stretch = 24 / 9

# lots of spurious overlaps detected when too extreme shifts allowed
stretch = c(2, 1.5, 1) # the stretch which lines up the start and end points.
min_num_overlapping_points = 4 # will only allow shifts which leave this many overlapping points after applying the stretch
shift_extreme  = 4 # the absolute maximum value which can be applied as a shift. Noticed that in the shuffled genes,
transformed.timecourse = 'Col0' # the name of the timecourse to apply the registratsion to (one of the names in the mean_df$accession column)
                                # which is loaded at line 133.

num.shuffled <- 1 #25 # for real, ran 40 jobs to get 40 * 25 random shuffled pairs for comparison.





print('********************')
print(paste0('output will be in ', outdir.string, '!'))
print('********************')


# setup flags for rescaling options
if (do.initial.rescale=='rescale') {
  initial.rescale <- TRUE
} else {
  initial.rescale <- FALSE
}
if (do.register.rescale=='rescale') {
  should.rescale <- TRUE
} else {
  should.rescale <- FALSE
}

#real.and.shuffled <- c('real', 1)

if (initial.rescale==TRUE) {
  print('********************')
  print('will rescale the data prior to registering, and register using this rescaled mean data!')
  print('********************')
}
if (should.rescale==TRUE){
  print('********************')
  print('will rescale the data when deciding optimal registration!')
  print('********************')
}


# directories to save graphs to
real.data.graph.dir <- paste0('../graphs/gene_registration/', outdir.string, '_real_data/')
shuffled.data.graph.dir <- paste0('../graphs/gene_registration/', outdir.string, '_shuffled_data/job_', jobNum, '/')
# directories to save real and shuffled expression data to
real.expression.dir <- paste0('../intermediate_data/gene_registration/', outdir.string, '_real_data/gene_expression/job_', jobNum, '/')
shuffled.expression.dir <- paste0('../intermediate_data/gene_registration/', outdir.string, '_shuffled_data/gene_expression/')
# directories to save real and shuffled distance data to
real.distance.dir <- paste0('../intermediate_data/gene_registration/', outdir.string, '_real_distance/job_', jobNum, '/')
shuffled.distance.dir <- paste0('../intermediate_data/gene_registration/', outdir.string, '_shuffled_distance/')


# somewhere to store the data.tables and graphs
if (!(dir.exists(real.expression.dir))) {
    dir.create(real.expression.dir, recursive=T)
    dir.create(shuffled.expression.dir, recursive=T)

    dir.create(real.distance.dir, recursive=T)
    dir.create(shuffled.distance.dir, recursive=T)

    dir.create(real.data.graph.dir, recursive=T)
    dir.create(shuffled.data.graph.dir, recursive=T)
}


## GET THE RAW DATA. MEAN EXPRESSION OF BIOLOGICAL REPS. TREAT BRASSICA GENES INDIVIDUALLY (DON'T SUM THEM).
# load the data expression data. Consider the brassica genes individually (don't sum)
L <- load_mean_df() # sumBrassica copy flag is within load_mean_df()
mean_df <- L[[1]] # will compare Col0 vs Ro18 based on accession column
all.data.df <- L[[2]]


#cut down to few genes for testing
test.genes <- unique(mean_df$locus_name)[1:200]
mean_df <- mean_df[mean_df$locus_name %in% test.genes,]
all.data.df <- all.data.df[all.data.df$locus_name %in% test.genes, ]


unshuffled.data <- copy(mean_df)
unshuffled.all.data <- copy(all.data.df)

real.and.shuffled <- c('real', 1:num.shuffled) # if just want to do for real=c('real'), or for random shuffles as well
i <- 'real'
for (i in real.and.shuffled) {

  print('*****************')
  print(paste0('on loop : ', i))
  print('*****************')


  if (i == 'real') {
    mean_df <- copy(unshuffled.data)
    all.data.df <- copy(unshuffled.all.data)
  } else {
    if (shuffle.type=='shuffle.expression') {
      print('shuffling expression timecourse for each gene')
      L <- shuffle_ro18_timepoints(mean_df, all.data.df)
      mean_df <- L[[1]]
      all.data.df <- L[[2]]
    } else if (shuffle.type=='shuffle.genes') {
      print('shuffling gene names')
      L <- shuffle_ro18_gene_names(mean_df, all.data.df)
      mean_df <- L[[1]]
      all.data.df <- L[[2]]
    } else {
      print(junk.because.wrong.shuffle.type)
    }
  }


  ## PREPARE, AND REGISTER AND SCALE THE DATA
  O <- prepare_scaled_and_registered_data(mean_df, all.data.df, stretch=stretch, initial.rescale, should.rescale, min_num_overlapping_points,
                                          shift_extreme, transformed.timecourse)

  mean_df <- O[['mean_df']] # mean_df is unchanged
  mean_df.sc <- O[['mean_df.sc']] # mean_df.sc : data is scaled(center=T, scale=T)
  imputed.mean_df <- O[['imputed.mean_df']] # imputed.mean_df is registered data, Col0 values imputed to make a single timepoint.
  all_shifts <- O[['all_shifts']] # all_shifts is data.table of score for each shift for each gene.
  model.comparison <- O[['model.comparison']]

  # sanity test plot
  # ggplot(imputed.mean_df[imputed.mean_df$locus_name=='BRAA02G015410.3C',],
  #        aes(x=shifted_time, y=mean.cpm, color=accession))+
  #   geom_point()+
  #   geom_line()

  #### CALCULATE THE DISTANCES BETWEEN THE SAMPLES ####
  O <- calculate_between_sample_distance(mean_df, mean_df.sc, imputed.mean_df)
  D.mean <- O[['D.mean']]
  D.scaled <- O[['D.scaled']]
  D.registered <- O[['D.registered']]
  D.scaled.onlyNR <- O[['D.scaled.onlyNR']]
  D.scaled.onlyR <- O[['D.scaled.onlyR']]
  D.registered.onlyR <- O[['D.registered.onlyR']]


  # Save the generated tables
  if (i=='real') {
    # save the expression info
    saveRDS(mean_df, file=paste0(real.expression.dir, 'mean_df.rds'))
    saveRDS(mean_df.sc, file=paste0(real.expression.dir, 'mean_df.sc.rds'))
    saveRDS(imputed.mean_df, file=paste0(real.expression.dir, 'imputed.mean_df.rds'))
    saveRDS(all_shifts, file=paste0(real.expression.dir, 'all_shifts.rds'))
    saveRDS(model.comparison, file=paste0(real.expression.dir, 'model.comparison.rds'))

    # save the distances calculated
    saveRDS(D.mean, file=paste0(real.distance.dir, 'D.mean.rds'))
    saveRDS(D.scaled, file=paste0(real.distance.dir, 'D.scaled.rds'))
    saveRDS(D.registered, file=paste0(real.distance.dir, 'D.registered.rds'))
    saveRDS(D.scaled.onlyNR, file=paste0(real.distance.dir, 'D.scaled.onlyNR.rds'))
    saveRDS(D.scaled.onlyR, file=paste0(real.distance.dir, 'D.scaled.onlyR.rds'))
    saveRDS(D.registered.onlyR, file=paste0(real.distance.dir, 'D.registered.onlyR.rds'))

  } else { # if the data timecourse has been shuffled
    # save the expression info
    saveRDS(mean_df, file=paste0(shuffled.expression.dir, 'mean_df_', jobNum, '_', i, '.rds'))
    saveRDS(mean_df.sc, file=paste0(shuffled.expression.dir, 'mean_df.sc_', jobNum, '_', i, '.rds'))
    saveRDS(imputed.mean_df, file=paste0(shuffled.expression.dir, 'imputed.mean_df_', jobNum, '_', i, '.rds'))
    saveRDS(all_shifts, file=paste0(shuffled.expression.dir, 'all_shifts_', jobNum, '_', i, '.rds'))
    saveRDS(model.comparison, file=paste0(shuffled.expression.dir, 'model.comparison_', jobNum, '_', i,'.rds'))


    # save the distances calculated
    saveRDS(D.mean, file=paste0(shuffled.distance.dir, 'D.mean_', jobNum, '_', i, '.rds'))
    saveRDS(D.scaled, file=paste0(shuffled.distance.dir, 'D.scaled_', jobNum, '_', i, '.rds'))
    saveRDS(D.registered, file=paste0(shuffled.distance.dir, 'D.registered_', jobNum, '_', i, '.rds'))
    saveRDS(D.scaled.onlyNR, file=paste0(shuffled.distance.dir, 'D.scaled.onlyNR', jobNum, '_', i,'.rds'))
    saveRDS(D.scaled.onlyR, file=paste0(shuffled.distance.dir, 'D.scaled.onlyR', jobNum, '_', i, '.rds'))
    saveRDS(D.registered.onlyR, file=paste0(shuffled.distance.dir, 'D.registered.onlyR', jobNum, '_', i, '.rds'))
  }


  # D.all <- do.call('rbind', list(D.mean, D.scaled, D.registered))
  # D.all$title <- factor(D.all$title, levels=c('mean expression', 'scaled mean expression',
  #                                             'registered & scaled mean expression'))
  # p.all <- make_heatmap_all(D.all, '')
  # ggsave(filename='../graphs/registration_plots/distance_plots_real_data_all.pdf',
  #        p.all,
  #        width=3,
  #        height=10)

  # replace self comparison with NA for scaled
  D.scaled$distance[grepl('Col0', D.scaled$x.sample) &
                          grepl('Col0', D.scaled$y.sample) ] <- NA
  D.scaled$distance[grepl('Ro18', D.scaled$x.sample) &
                          grepl('Ro18', D.scaled$y.sample) ] <- NA

  D.scaled.onlyNR$distance[grepl('Col0', D.scaled.onlyNR$x.sample) &
                      grepl('Col0', D.scaled.onlyNR$y.sample) ] <- NA
  D.scaled.onlyNR$distance[grepl('Ro18', D.scaled.onlyNR$x.sample) &
                      grepl('Ro18', D.scaled.onlyNR$y.sample) ] <- NA

  D.scaled.onlyR$distance[grepl('Col0', D.scaled.onlyR$x.sample) &
                      grepl('Col0', D.scaled.onlyR$y.sample) ] <- NA
  D.scaled.onlyR$distance[grepl('Ro18', D.scaled.onlyR$x.sample) &
                      grepl('Ro18', D.scaled.onlyR$y.sample) ] <- NA



  # replace self comparison with NA for registered
  D.registered$distance[grepl('Col0', D.registered$x.sample) &
                        grepl('Col0', D.registered$y.sample) ] <- NA
  D.registered$distance[grepl('Ro18', D.registered$x.sample) &
                          grepl('Ro18', D.registered$y.sample) ] <- NA

  D.registered.onlyR$distance[grepl('Col0', D.registered.onlyR$x.sample) &
                          grepl('Col0', D.registered.onlyR$y.sample) ] <- NA
  D.registered.onlyR$distance[grepl('Ro18', D.registered.onlyR$x.sample) &
                          grepl('Ro18', D.registered.onlyR$y.sample) ] <- NA

  p.all <- make_data_heatmaps(D.mean, D.scaled, D.registered, D.scaled.onlyNR, D.scaled.onlyR, D.registered.onlyR)
  if (i=='real') {
    ggsave(filename=paste0(real.data.graph.dir, 'job_', jobNum, '_distance_plots_real_data_all_sepScale.pdf'),
           p.all,
           width=12,
           height=16)
  } else {
    ggsave(filename=paste0(shuffled.data.graph.dir, i, '_distance_plots_shuffled_data_all_sepScale.pdf'),
           p.all,
           width=12,
           height=16)
  }
}

