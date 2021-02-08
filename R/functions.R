#file.type <- 'mean.df'

#' @export
load_shuffled_data <- function(shuffled.data.dir, file.type) {
  # file.type is type of file loading:
  # "comparison": for model.comparison files
  # "mean.sc": for mean.df.sc files
  # "mean.df" : for mean.df files
  # "imputed.mean.df" for them
  # "shifts" for all.shifts

  allowed.types <- c('comparison', 'mean.sc', 'mean.df', 'imputed.mean.df', 'shifts')
  if (!(file.type %in% allowed.types)) {
    print(paste0('file.type must be one of : ', paste0(allowed.types, collapse=', ')))
    stop()
  }
  # format file.type to searchable pattern - handle "."s
  if (file.type=='mean.sc') {
    file.type <- 'mean\\.df\\.sc'
  } else if (file.type=='mean.df') {
    file.type <-  '^mean\\.df'
  } else if (file.type=='imputed.mean.df') {
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

#data.dir <- shuffled.data.dir
#' @export
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

#all.rep.data <- all.data.df
#' @export
scale_all_rep_data <- function(mean.df, all.rep.data, scale.func) {
  # apply the same scaling which done to the mean expression data
  # to all the reps.
  # (subtract mean, and divide by sd), using the values for the mean data
  # as this is what was used to find the best shift.

  # calculate the summary stats to use for the rescaling
  gene.expression.stats <- unique(mean.df[,
                                          .(mean_val=mean(mean.cpm),
                                            sd_val=stats::sd(mean.cpm)),
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

  # ggplot2::ggplot(mean.df[mean.df$locus_name=='BRAA01G000040.3C', ], ) +
  #  ggplot2::aes(x=timepoint, y=mean.cpm, color=accession
  #  ggplot2::geom_point()


  return(out)
}

#' @export
my.scale <- function(v) {
  return(v / max(v))
}

#' @export
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

#' @export
load_mean.df <- function() {

  #setwd('/Volumes/Research-Projects/bravo/alex/BRAVO_rna-seq/scripts/')
  rds_file <- 'ro18_chiifu_apex' # don't include the .rds # the name of the brassica data to load
  sumBrassicas <- F # if false use seperate brassica genes, and compare to repeated arabidopsis genes. If true, sume copies of each brassica and compare to arabidopsis
  #datapath <- paste0('../final_data/rds/', rds_file, '.rds')


  #### specify the genes to be used in the comparison: ####
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

# calculate the comparison stats
#' @export
calc.AIC <- function(logL, num.params) {
  return((-2*logL) + 2*num.params)
}

#' @export
calc.BIC <- function(logL, num.params, num.obs) {
  return((-2*logL) + log(num.obs) * num.params)
}


#' @export
change.accession.names <- function(mean.df, all.data.df, transformed.timecourse) {
  # set the "transformed.timecourse" accession to "Col0", and the other one to "Ro18"

  # error checking
  if (length(unique(mean.df$accession)) != 2) {
    stop('Error in change.accession.names() : comparison must be made between two accessions!')
  }

  # store these, to rename at the end
  original.transformed.timecourse.name <- transformed.timecourse
  original.other.accession.name <- as.character(unique(mean.df$accession[mean.df$accession!=transformed.timecourse]))

  # change mean.df
  new.mean.df.accession <- mean.df$accession
  new.mean.df.accession[mean.df$accession==transformed.timecourse] <- 'Col0'
  new.mean.df.accession[mean.df$accession!=transformed.timecourse] <- 'Ro18'
  mean.df$accession <- new.mean.df.accession

  # change all.data.df
  new.all.data.df.accession <- all.data.df$accession
  new.all.data.df.accession[all.data.df$accession==transformed.timecourse] <- 'Col0'
  new.all.data.df.accession[all.data.df$accession!=transformed.timecourse] <- 'Ro18'
  all.data.df$accession <- new.all.data.df.accession

  return(list('mean.df'=mean.df,
              'all.data.df'=all.data.df,
              'original.transformed.accession.name'=original.transformed.timecourse.name,
              'original.other.accession.name'=original.other.accession.name))
}


# stretches=c(1, 1.5, 2.0)
# initial.rescale=FALSE
# do.rescale=FALSE
# min.num.overlapping.points = 4
# test.genes <- unique(mean.df$locus_name)[1:101]
# #test.genes <- 'MSTRG.10244'
# mean.df <- mean.df[mean.df$locus_name %in% test.genes,]
# all.data.df <- all.data.df[all.data.df$locus_name %in% test.genes,]
# shift.extreme=4
# transformed.timecourse <- 'Ro18'

#' @export
prepare_scaled_and_registered_data <- function(mean.df, all.data.df, stretches,
                                               initial.rescale, do.rescale,
                                               min.num.overlapping.points, shift.extreme,
                                               transformed.timecourse) {
  ## APPLY NORMALISATION OF EXPRESSION FOR EACH GENE ACROSS ALL TIMEPOINTS ##

  # hardcoded all the functions to use 'Col0', and 'Ro18'. Rather than fix all instances
  # of that, just temporarily rename them here, and turn back at the end.
  # will apply stretch, and shift to the "transformed.timecourse" accession (which should be the quicker one...)
  L <- change.accession.names(mean.df, all.data.df, transformed.timecourse)
  mean.df <- L[['mean.df']]
  all.data.df <- L[['all.data.df']]
  original.transformed.accession <- L[['original.transformed.accession.name']]
  original.other.accession <- L[['original.other.accession.name']]



  mean.df.sc <- data.table::copy(mean.df)
  # specify what kind of scaling
  mean.df.sc[, sc.mean.cpm:=scale(mean.cpm, scale=TRUE, center=TRUE), by=.(locus_name, accession)]
  #mean.df.sc[, sc.mean.cpm:=my.scale(mean.cpm), by=.(locus_name, accession)]

  ## APPLY INDIVIDUAL SHIFTING ######
  # optimise transformations applied to arabidopsis gene profile to map onto the brassicas - shift in x direction, using mean for mapping,
  # and only rescale expression using the candidate timepoints in common.

  # specify which data to use for registering.
  # whether prior rescaled mean, or mean data should be used for registration
  if (initial.rescale==TRUE) {
    # apply rescale to mean.df prior to registration
    to.shift.df <- data.table::copy(mean.df.sc)
    to.shift.df$mean.cpm <- to.shift.df$sc.mean.cpm
    to.shift.df$sc.mean.cpm <- NULL

    # apply THE SAME rescale to all.data.df prior to registration
    #all.data.df <- scale_all_rep_data(mean.df, all.data.df, 'my.scale')
    all.data.df <- scale_all_rep_data(mean.df, all.data.df, 'scale')


    # sanity plot that rescale all data worked
    # ggplot2::ggplot(all.data.df[all.data.df$locus_name=='BRAA01G000040.3C'])+
    #   ggplot2::aes(x=timepoint, y=mean.cpm, color=accession)
    #   ggplot2::geom_point()
  } else {
    to.shift.df <- data.table::copy(mean.df)
  }
  # ggplot2::ggplot(to.shift.df[to.shift.df$locus_name=='BRAA03G004600.3C'])+
  #   ggplot2::aes(x=timepoint, y=mean.cpm, color=accession)
  #   ggplot2::geom_point()
  # tst <- all.data.df
  # ggplot2::ggplot(tst[tst$locus_name=='BRAA01G000040.3C'])+
  #   ggplot2::aes(x=timepoint, y=mean.cpm, color=accession) +
  #   ggplot2::geom_point()

  # calculate the best registration. Returns all tried registrations, best stretch and shift combo,
  # and AIC/BIC stats for comparison of best registration model to seperate models for expression of
  # each gene in Ro18 and Col0.
  L <- get_best_stretch_and_shift(to.shift.df, all.data.df, stretches, do.rescale, min.num.overlapping.points, shift.extreme)
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
  shifted.mean.df <- apply_shift_to_registered_genes_only(to.shift.df, best_shifts, model.comparison.dt)
  # shifted.mean.df <- apply_best_shift(to.shift.df, best_shifts) # can be NA if exactly tied for what the best shift was

  # GOI <- 'MSTRG.11237'
  # ggplot2::ggplot(all.data.df[all.data.df$locus_name==GOI])+
  #   ggplot2::aes(x=timepoint, y= mean.cpm, color=accession) +
  #   ggplot2::geom_point()
  #
  # #sanity plot that done right
  # ggplot2::ggplot(shifted.mean.df[shifted.mean.df$locus_name==GOI])+
  #   ggplot2::aes(x=shifted.time, y=mean.cpm, color=accession) +
  #   ggplot2::geom_point()+
  #   ggplot2::geom_line()


  # impute arabidopsis values at times == to the observed brassica points for each shifted arabidopsis gene
  # so can compare using heatmap.
  # arabidopsis curves are the ones that been shifted around. Linear impute values for these
  # curves so that brassica samples can be compared to an arabidopsis point.
  imputed.mean.df <- impute_arabidopsis_values(shifted.mean.df)

  #sanity plot that done right
  # ggplot2::ggplot(shifted.mean.df[shifted.mean.df$locus_name=='BRAA01G001540.3C'])+
  #   ggplot2::aes(x=shifted.time, y=mean.cpm, color=accession) +
  #   ggplot2::geom_point()+
  #   ggplot2::geom_line()
  # ggplot2::ggplot(imputed.mean.df[imputed.mean.df$locus_name=='BRAA01G001540.3C'])+
  #   ggplot2::aes(x=shifted.time, y=mean.cpm, color=accession)+
  #   ggplot2::geom_point()+
  #   ggplot2::geom_line()


  # fix the accession names to the ones actually passed in:
  mean.df <- fix.accessions(mean.df, original.transformed.accession, original.other.accession)
  mean.df.sc <- fix.accessions(mean.df.sc, original.transformed.accession, original.other.accession)
  imputed.mean.df <- fix.accessions(imputed.mean.df, original.transformed.accession, original.other.accession)


  OUT <- list('mean.df'=mean.df,
              'mean.df.sc'=mean.df.sc,
              'imputed.mean.df'=imputed.mean.df,
              'all.shifts'=all_shifts,
              'model.comparison'=model.comparison.dt)
}

#' @export
fix.accessions <- function(df, original.transformed.accession, original.other.accession) {
  # swap Col0 with original.transformed.accession, and Ro18 with original.other.accession
  new.df.accession <- df$accession
  new.df.accession[df$accession=='Col0'] <- original.transformed.accession
  new.df.accession[df$accession=='Ro18'] <- original.other.accession
  df$accession <- new.df.accession

  return(df)
}


#' @export
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
      min.shift <- unique(cand.shifts[cand.shifts==min(cand.shifts)])
      is.best[df$shift != min.shift] <- FALSE
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
# do.rescale
# min.num.overlapping.points
# shift.extreme
#' @export
get_best_stretch_and_shift <- function(to.shift.df, all.data.df, stretches, do.rescale, min.num.overlapping.points, shift.extreme) {
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
    # if do.rescale=T, is rescaled by the mean FOR THE OVERLAPPING POINTS. (but not by the SD.)

    # ggplot2::ggplot(to.shift.df[to.shift.df$locus_name=='MSTRG.12467',])+
    #   ggplot2::aes(x=timepoint, y=mean.cpm, color=accession) +
    #   ggplot2::geom_point()

    all_shifts <- calculate_all_best_shifts(to.shift.df, stretch_factor=stretch, do.rescale, min.num.overlapping.points, shift.extreme)

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


#' @export
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
  seperate.dt$shifted.time <- seperate.dt$stretched.time.delta + 11 # add eleven, as this is done for the registered genes
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

#' @export
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
  # ggplot2::ggplot(tst[tst$locus_name=='BRAA01G000040.3C'])+
  #   ggplot2::aes(x=shifted.time, y=mean.cpm, color=accession) +
  #   ggplot2::geom_point()
  # ggplot2::ggplot(tst[tst$locus_name=='BRAA01G000040.3C'])+
  #   ggplot2::aes(x=timepoint, y=mean.cpm, color=accession) +
  #   ggplot2::geom_point()

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

  out <- data.table::data.table(data.frame('gene'=genes, 'seperate.AIC'=out.sepAIC, 'registered.AIC'=out.combAIC,
                                           'seperate.BIC'=out.sepBIC, 'registered.BIC'=out.combBIC))
  return(out)
}

# curr.sym <- 'BRAA01G001320.3C'
# all.data.df <- all.data.df
# is.testing <- TRUE
#' @export
compare_registered_to_unregistered_model <- function(curr.sym, all.data.df, is.testing) {
  # compare the overlapping timepoints in brassica and arabidopsis after the best registration,
  # and without registration (use the same timepoints for both models).
  # use the stretched data for both models, whether considering as registered or not.


  curr.data.df <- all.data.df[all.data.df$locus_name==curr.sym]

  # print('line 662')
  # print(curr.data.df)

  # flag the timepoints to be used in the modelling, only the ones which overlap!
  curr.data.df <- get_compared_timepoints(curr.data.df)

  # ggplot2::ggplot(curr.data.df)+
  #   ggplot2::aes(x=shifted.time, y=mean.cpm, shape=is.compared, color=accession)+
  #   ggplot2::geom_point()

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


  ara.fit <- stats::lm(mean.cpm~splines::bs(shifted.time, df=num.spline.params, degree=3), data=ara.spline.data)
  bra.fit <- stats::lm(mean.cpm~splines::bs(shifted.time, df=num.spline.params, degree=3), data=bra.spline.data)
  combined.fit <- stats::lm(mean.cpm~splines::bs(shifted.time, df=num.spline.params, degree=3), data=combined.spline.data)
  # calculate the log likelihoods
  ara.logLik <- stats::logLik(ara.fit)
  bra.logLik <- stats::logLik(bra.fit)
  seperate.logLik <- ara.logLik + bra.logLik # logLikelihoods, so sum
  combined.logLik <- stats::logLik(combined.fit)

  # calculate the comparison.stats - - AIC, BIC, smaller is better!
  # 2*num.spline.params as fitting seperate models for Ara * Col
  seperate.AIC <- calc.AIC(seperate.logLik, 2*num.spline.params)
  combined.AIC <- calc.AIC(combined.logLik, num.spline.params+num.registration.params)

  seperate.BIC <- calc.BIC(seperate.logLik, 2*num.spline.params, num.obs)
  combined.BIC <- calc.BIC(combined.logLik, num.spline.params+num.registration.params, num.obs)


  if (is.testing==TRUE) {
    ara.pred <- stats::predict(ara.fit)
    ara.pred.df <- unique(data.frame('shifted.time'=ara.spline.data$shifted.time,
                                     'mean.cpm'=ara.pred, 'accession'='Col0'))
    bra.pred <- stats::predict(bra.fit)
    bra.pred.df <- unique(data.frame('shifted.time'=bra.spline.data$shifted.time,
                                     'mean.cpm'=bra.pred, 'accession'='Ro18'))

    combined.pred <- stats::predict(combined.fit)
    combined.pred.df <- unique(data.frame('shifted.time'=combined.spline.data$shifted.time,
                                          'mean.cpm'=combined.pred, 'accession'='registered'))
    spline.df <- rbind(ara.pred.df, bra.pred.df, combined.pred.df)

    ggplot2::ggplot(data=combined.spline.data)+
      ggplot2::aes(x=shifted.time, y=mean.cpm, colour=accession) +
      ggplot2::geom_point()+
      ggplot2::geom_line(data=spline.df)+
      ggplot2::ggtitle(paste0(curr.sym, ' : sep AIC:combo AIC=', round(seperate.AIC), ':', round(combined.AIC),
                              ', sep BIC: combo BIC=', round(seperate.BIC), ':', round(combined.BIC)))
    ggplot2::ggsave(paste0('./testing/fitted_splines/', curr.sym, '_', max(ara.pred.df$shifted.time), '.pdf'))
  }

  return(list(seperate.AIC, combined.AIC, seperate.BIC,combined.BIC))
}




#' @export
make_heatmap <- function(D, ylabel, y.axis.fontsize=6) {
  D$x.sample <- factor(D$x.sample, levels=unique(sort(D$x.sample)))
  D$y.sample <- factor(D$y.sample, levels=unique(sort(D$y.sample)))

  p <- ggplot2::ggplot(D)+
    ggplot2::aes(x=x.sample, y=y.sample, fill=log(distance)) +
    ggplot2::geom_tile()+
    #scale_fill_viridis()+
    #theme_classic()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90,
                                                       size=6),
                   axis.text.y = ggplot2::element_text(size=y.axis.fontsize),
                   plot.title = ggplot2::element_text(hjust=0.5, size=10),
                   plot.margin = ggplot2::margin(0,0,-10,0),
                   panel.background = ggplot2::element_blank(),
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
    ggplot2::facet_wrap(~title, nrow=1)+
    ggplot2::ylab(ylabel)+
    ggplot2::xlab('')
  # ggplot2::ggtitle(title)

  return(p)
}

# mean.df <- real.mean.df
# mean.df.sc <- real.sc.df
# imputed.mean.df <- imputed.mean.df
#' @export
calculate_between_sample_distance <- function(mean.df, mean.df.sc, imputed.mean.df) {

  ### convert all to wide format ready for distance calculation

  # mean.df
  sample.id.cols <- c('accession','timepoint')
  gene.col <- c('locus_name')
  expression.col <- 'mean.cpm'
  mean.dt.w <- reformat_for_distance_calculation(mean.df, sample.id.cols, gene.col, expression.col)

  # normalised mean.df
  sample.id.cols <- c('accession','timepoint')
  gene.col <- c('locus_name')
  expression.col <- c('sc.mean.cpm')
  mean.dt.sc.w <- reformat_for_distance_calculation(mean.df.sc, sample.id.cols, gene.col, expression.col)

  # imputed.mean.df - all genes
  sample.id.cols <- c('accession','shifted.time')
  gene.col <- c('locus_name')
  expression.col <- c('mean.cpm')
  imputed.mean.dt.w <- reformat_for_distance_calculation(imputed.mean.df, sample.id.cols, gene.col, expression.col)

  # same, but for subsets of REGISTERED / NOT REGISTERED genes.
  # distance between samples, only using genes which are found best model is not registered
  not.registered.genes <- unique(imputed.mean.df$locus_name[imputed.mean.df$is.registered==FALSE])
  mean.dt.sc.w.not.registered <- mean.dt.sc.w[mean.dt.sc.w$locus_name %in% not.registered.genes,]
  # distance between samples, only using genes which are found best when ARE registered
  registered.genes <- unique(imputed.mean.df$locus_name[imputed.mean.df$is.registered==TRUE])
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


# unrandom.mean.df <- mean.df
# unrandom.all.df <- all.data.df
#' @export
shuffle_ro18_timepoints <- function(unrandom.mean.df, unrandom.all.df) {
  # shuffle the timepoints for each ro18 gene

  # split the mean.df
  col.df <- unrandom.mean.df[unrandom.mean.df$accession=='Col0',]
  ro18.df <- unrandom.mean.df[unrandom.mean.df$accession=='Ro18',]
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

  mean.df <- rbind(col.df, shuffled.ro18.df)
  all.df <- rbind(col.all.df, shuffled.ro18.all.df)
  return(list(mean.df, all.df))
}

#' @export
shuffle_ro18_gene_names <- function(mean.df, out.all.df) {
  # shuffle the identities of the genes in the brassica
  # can't just do shuffle, becuase need to preserve which timepoints are from the same gene
  out.mean.df <- data.table::copy(mean.df)
  out.all.df <- data.table::copy(out.all.df)

  # make the gene lookup table for the same shuffled genes for both
  brassica.genes <- unique(out.mean.df$locus_name[out.mean.df$accession=='Ro18'])
  shuffled.genes <- sample(brassica.genes)
  shuffle.gene.lookup <- data.table::data.table(data.frame('gene.id'=brassica.genes, 'shuffled.id'=shuffled.genes))

  # change the gene names for the mean.df
  out.mean.df <- swap_gene_names(out.mean.df, shuffle.gene.lookup)
  # change the gene names for the all.df
  out.all.df <- swap_gene_names(out.all.df, shuffle.gene.lookup)

  return(list(out.mean.df, out.all.df))
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


#' @export
make_data_heatmaps <- function(D.mean, D.scaled, D.registered, D.scaled.NR, D.scaled.R, D.registered.R) {

  p.mean <- make_heatmap(D.mean, 'mean expression')
  p.scaled <- make_heatmap(D.scaled, 'scaled mean expression')
  p.registered <- make_heatmap(D.registered, 'registered & scaled mean expression')

  p.scaled.NR <- make_heatmap(D.scaled.onlyNR, 'scaled mean expression (only not registered genes)')
  p.scaled.R <- make_heatmap(D.scaled.onlyR, 'scaled mean expression (only registered genes)')
  p.registered.R <- make_heatmap(D.registered.R, 'registered & scaled mean expression (only  registered genes)')


  p.all <- cowplot::plot_grid(p.mean, p.scaled.NR, p.scaled, p.scaled.R, p.registered, p.registered.R, ncol=2)

  return(p.all)
}

#' @export
make_heatmap_all <- function(D, title) {
  D$x.sample <- factor(D$x.sample, levels=unique(sort(D$x.sample)))
  D$y.sample <- factor(D$y.sample, levels=unique(sort(D$y.sample)))

  p <- ggplot2::ggplot(D)+
    ggplot2::aes(x=x.sample, y=y.sample, fill=log(distance)) +
    ggplot2::geom_tile()+
    #viridis::scale_fill_viridis()+
    ggplot2::theme_classic()+
    ggplot2::facet_wrap(~title, ncol=1, scales='free')+
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
    ggplot2::ylab('')+
    ggplot2::xlab('')+
    ggplot2::ggtitle(title)

  return(p)
}

#' @export
calculate_pairwise_sample_distance_simple <- function(dt) {
  # dt is the wide format expression of the two samples to be compared.
  # only genes which have data in both samples are considered (only
  # relevant for registration case).

  # filter to genes with data in both.
  dt <- stats::na.omit(dt)
  # calculate distance
  dt$sq.diff <- (dt[, 1] - dt[, 2])^2

  d <- mean(dt$sq.diff)
  return(d)
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

#brassica_name <- 'ro18_chiifu_apex'
#' @export
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
#models_o_I <- 'filt_models'

#' @export
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

# sample.id.cols <- c('accession','delta.time')
# gene.col <- c('locus_name')
# expression.col <- 'mean.cpm'
# dt <- mean.df
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
#' @export
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

#' @export
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



#' @export
impute_arabidopsis_values <- function(shifted.mean.df) {
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
  # ggplot2::ggplot(shifted.mean.df[shifted.mean.df$locus_name=='BRAA01G001090.3C', ])+
  #   ggplot2::aes(x=shifted.time, y=mean.cpm, color=accession)+
  #   ggplot2::geom_point()
  # # an unregistered one
  # ggplot2::ggplot(shifted.mean.df[shifted.mean.df$locus_name=='BRAA01G003140.3C', ])+
  #   ggplot2::aes(x=shifted.time, y=mean.cpm, color=accession)+
  #   ggplot2::geom_point()

  # The imputed col0 times going to extimate gene expression for
  imputed.timepoints <- round(seq(min(shifted.mean.df$shifted.time), max(shifted.mean.df$shifted.time)))

  out.list <- list()
  out.list <- c(out.list, list(shifted.mean.df[shifted.mean.df$accession=='Ro18']))

  curr.gene <- 'BRAA01G001090.3C' #unique(shifted.mean.df$locus_name)[17]
  count <- 0
  for (curr.gene in unique(shifted.mean.df$locus_name)) {
    if (count %% 100 == 0) {
      print(paste0(count, ' / ', length(unique(shifted.mean.df$locus_name))))
    }

    # get the current gene expression data
    curr.df <- shifted.mean.df[shifted.mean.df$locus_name==curr.gene, ]

    # skip over this one, if not assigned a best shift value, because brassica gene not expressed
    # if (all(is.na(curr.df$mean.cpm))) {
    #   next
    # }

    ara.df <- curr.df[curr.df$accession=='Col0',]
    bra.df <- curr.df[curr.df$accession=='Ro18',]

    interp.ara.df <- data.table::data.table(data.frame('locus_name'=curr.gene, 'accession'='Col0', 'tissue'='apex', 'timepoint'=NA,
                                                       'stretched.time.delta'= NA, 'shifted.time'=imputed.timepoints,
                                                       'is.registered'= unique(ara.df$is.registered)[1]))

    # for each brassica timepoint, interpolate the comparible arabidopsis expression
    # by linear interpolation between the neighbouring 2 ara values. If not between 2 ara values
    # because shifted outside comparible range, set to NA
    interp.ara.df$mean.cpm <- sapply(imputed.timepoints, interpolate_brassica_comparison_expression, bra.dt=ara.df)


    # sanity testing - line is interpolated.
    # ggplot2::ggplot(curr.df)+
    #   ggplot2::aes(x=shifted.time, y=mean.cpm, color=accession)+
    #   ggplot2::geom_point()+
    #   ggplot2::geom_line(data=interp.ara.df)

    out.list <- c(out.list, list(interp.ara.df))
    count <- count+1
  }
  out.df <- do.call('rbind', out.list)
  return(out.df)
}

# test <- mean.df
# stretch_factor
# min.num.overlapping.points
# shift.extreme
#' @export
get_extreme_shifts_for_all <- function(test, stretch_factor, min.num.overlapping.points, shift.extreme) {
  # wrapper for calc_extreme_shifts to be able to move it out of the loop so don't calculate for every gene.

  #min.num.overlapping.points <- 5 # bound the extreme allowed shifts, such than at least this many timepoints are being compared for both accessions.
  # cut data.table to a single gene
  curr_sym <- unique(test$locus_name)[1]
  test <- test[test$locus_name==curr_sym, ]

  # transform timepoint to be time from first timepoint
  test[, delta.time:=timepoint - min(timepoint), by=.(accession)]
  # apply stretch_factor to the arabidopsis, leave the rapa as is
  test$delta.time[test$accession=='Col0'] <- test$delta.time[test$accession=='Col0']*stretch_factor

  # calculate min shift and max time shift, which still allows overlap of at least 5 times to be compared from whichever accession will be considering fewer timepoints from.
  # Shift is applied to the arabidopsis - so the 5th largest arabidopsis time is the biggest -ve shift can be applied
  # and the biggest shift which can be applied is to make the 5th smallest arabidopsis time == largest brassica time
  #data.table::setorder(test, delta.time)
  #min.shift <- min(test$delta.time[test$accession=='Ro18']) - test$delta.time[test$accession=='Col0'][length(test$delta.time[test$accession=='Col0'])-4]
  #max.shift <- max(test$delta.time[test$accession=='Ro18']) - test$delta.time[test$accession=='Col0'][5]

  M <- calc_extreme_shifts(test, min.num.overlapping.points, shift.extreme)
  return(M)
}

# for all genes, in mean.df, get the scores, and the shifts which result in them after
# stretching arabidopsis by "stretch_factor, and applying shifts forward and
# backward, with the extremes defined by allowing 5 overlapping points for comparison.
# mean.df <- to.shift.df
# stretch_factor <- 2
# do.rescale=F
#' @export
calculate_all_best_shifts <- function(mean.df, stretch_factor, do.rescale, min.num.overlapping.points, shift.extreme) {
  symbols <- c()
  num_points <- c()
  #curr_sym <- 'TT16'
  all.scores.list <- rep(list(0), length(unique(mean.df$locus_name)))
  length(unique(mean.df$locus_name))


  # get the extreme shifts which can be applied to the genes
  M <- get_extreme_shifts_for_all(mean.df, stretch_factor, min.num.overlapping.points, shift.extreme)
  min.shift <- M[[1]]
  max.shift <- M[[2]]

  count <- 0
  curr_sym <- unique(mean.df$locus_name)[2]

  i = 1
  for (i in 1:length(unique(mean.df$locus_name))) {
    #for (curr_sym in unique(mean.df$locus_name)) {
    curr_sym <- unique(mean.df$locus_name)[i]
    if (count %% 100 == 0) {
      print(paste0(count, ' / ', length(unique(mean.df$locus_name))))
    }

    # out is mean SSD between arabidopsis, and interpolated brassica (interpolated between 2 nearest points)
    # ggplot2::ggplot(mean.df[mean.df$locus_name==curr_sym,])+
    #   ggplot2::aes(x=timepoint, y=mean.cpm, color=accession) +
    #   ggplot2::geom_point()

    ### get "score" for all the candidate shifts - score is mean error / brassica expression for compared points.
    ### if timepoints don't line up, brassica value is linearly imputed
    out <- get_best_shift_new(curr_sym, mean.df, stretch_factor, do.rescale, min.shift, max.shift, testing=FALSE)

    best_shift <- out$shift[out$score==min(out$score)]
    if (length(best_shift) > 1) {
      if (max(out$score)=='Inf') { # can get inf score if brassica gene note expressed in the comparison
        next
      } else {
        # if ties for the best shift applied, apply the smaller absolute one
        best_shift <- best_shift[abs(best_shift) == min(abs(best_shift))]
      }
    }

    all.scores <- out
    all.scores.list[[i]] <- all.scores
    symbols <- c(symbols, curr_sym)
    #num_points <- c(num_points, best_shift)
    #print(best_shift)

    count <- count + 1
  }
  #shift_results <- data.frame(symbol=symbols, num.points = num_points)
  all.scores.df <- do.call('rbind', all.scores.list)

  return(all.scores.df)
}

# wroking on new implementation of the funciton
# curr_sym <- 'BRAA02G015410.3C'
#test <- mean.df
# stretch_factor <- 1
# do.rescale <- TRUE
#testing=T

#' @export
get_best_shift_new <- function(curr_sym, test, stretch_factor, do.rescale, min.shift, max.shift, testing=FALSE) {
  # for the current gene, and current stretch_factor, calculate the score for all
  # shifts, and return the scores for all as a table, and the value of the optimal shift.

  # Shift extremes are defined s.t. at least 5 points are compared.

  # do.rescale == TRUE, means apply "scale" to compared points for each shift. ==FALSE, means use original mean expression data

  num.shifts <- 10 # the number of different shifts to be considered.


  test <- test[test$locus_name==curr_sym, ]

  # transform timepoint to be time from first timepoint
  test[, delta.time:=timepoint - min(timepoint), by=.(accession)]
  # apply stretch_factor to the arabidopsis, leave the rapa as is
  test$delta.time[test$accession=='Col0'] <- test$delta.time[test$accession=='Col0']*stretch_factor

  all.scores <- rep(0, num.shifts)
  all.ara.mean <- rep(0, num.shifts)
  all.bra.mean <- rep(0, num.shifts)
  all.ara.sd <- rep(0, num.shifts)
  all.bra.sd <- rep(0, num.shifts)

  all.shifts <- seq(min.shift, max.shift, length.out=num.shifts)
  if (!(0 %in% all.shifts)) {
    all.shifts <- c(all.shifts, 0) # include 0 shift in candidates.
  }
  i=1
  for (i in 1:length(all.shifts)) {

    curr.shift <- all.shifts[i]

    #print('line 1676')
    #print(curr.shift)

    # shift the arabidopsis expression timeings
    test$shifted.time <- test$delta.time
    test$shifted.time[test$accession=='Col0'] <- test$delta.time[test$accession=='Col0'] + curr.shift

    #### test plot - of shifted, UNNORMALISED gene expression
    # if (testing==TRUE) {
    #   p <- ggplot2::ggplot(test)+
    #     ggplot2::aes(x=shifted.time, y=mean.cpm, color=accession)+
    #     ggplot2::geom_point()+
    #     ggplot2::ggtitle(paste0('shift : ', curr.shift))
    #   p
    #   ggplot2::ggsave(paste0('./testing/', curr.shift, '.pdf'))
    # }


    # cut down to just the arabidopsis and brassica timepoints which compared
    test <- get_compared_timepoints(test)
    compared <- test[test$is.compared==TRUE, ]

    # renormalise expression using just these timepoints?
    if (do.rescale==TRUE) {
      # record the mean and sd of the compared points, used for rescaling
      # in "apply shift" function
      ara.mean <- mean(compared$mean.cpm[compared$accession=='Col0'])
      bra.mean <- mean(compared$mean.cpm[compared$accession=='Ro18'])
      ara.sd <- stats::sd(compared$mean.cpm[compared$accession=='Col0'])
      bra.sd<- stats::sd(compared$mean.cpm[compared$accession=='Ro18'])

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
      p <- ggplot2::ggplot(compared)+
        ggplot2::aes(x=shifted.time, y=mean.cpm, color=accession)+
        ggplot2::geom_point()+
        ggplot2::ggtitle(paste0('shift : ', curr.shift))
      ggplot2::ggsave(paste0('./testing/',stretch_factor, '-', curr.shift, '.pdf'))
    }

    # for each arabidopsis timepoint, linear interpolate between the two nearest brassica timepoints
    ara.compared <- compared[compared$accession=='Col0']
    bra.compared <- compared[compared$accession=='Ro18']

    ara.compared$pred.bra.expression <- sapply(ara.compared$shifted.time, interpolate_brassica_comparison_expression, bra.dt=bra.compared)

    # calculate the score, using the (interpolated) predicted.bra.expression, and the observed arabidopsis expression
    # score = mean ((observed - expected)**2 )
    score <- calc_score(ara.compared$mean.cpm, ara.compared$pred.bra.expression)

    if (is.na(score)) {
      print('error in get_best_shift_new(): got a score of NA for gene:')
      print(curr_sym)
      print(paste0('with curr.shift=', curr.shift))
      stop()
    }

    all.scores[i] <- score
    all.ara.mean[i] <- ara.mean
    all.bra.mean[i] <- bra.mean
    all.ara.sd[i] <- ara.sd
    all.bra.sd[i] <- bra.sd
  }

  out <- data.table::data.table(data.frame('gene'=curr_sym, 'stretch'=stretch_factor, 'shift'=all.shifts, 'score'=all.scores,
                                           'ara.compared.mean'=all.ara.mean, 'bra.compared.mean'=all.bra.mean,
                                           'ara.compared.sd'=all.ara.sd, 'bra.compared.sd'=all.bra.sd))
  return(out)
}


#' @export
calc_num_overlapping_points <- function(shift, original) {
  # calculate the number of overlapping points for the species with the fewer overlapping points if the current "shift" is
  # applied to the col0 delta timepoints.
  original$shifted.time[original$accession=='Col0'] <- original$delta.time[original$accession=='Col0'] + shift
  original <- get_compared_timepoints(original)
  original[, num.compared:=sum(is.compared), by=.(accession)]

  return(min(original$num.compared))
}

#' @export
calc_extreme_shifts <- function(test, min.num.overlapping.points, shift.extreme) {
  # calculate the minimum and maximum shifts can apply to Col-0 after the stretch transformation, whilst
  # preserving the criteria that at least min.num.overlapping.points are being compared from both accessions.

  original <- data.table::copy(test)
  original$shifted.time <- original$delta.time

  # print('line 1803')
  # print(original)

  # -ve extreme shift will be -1*exactly the difference between 1 of the stretched Col0 timepoints, and the smallest ro18 timepoint
  # +ve extreme will be the difference between 1 of the col0 timepoints, and the maximum Ro18 timepoint
  neg.extreme.candidates <- -1*(original$delta.time[original$accession=='Col0'] - min(original$delta.time[original$accession=='Ro18']))
  pos.extreme.candidates <- max(original$delta.time[original$accession=='Ro18']) - original$delta.time[original$accession=='Col0']

  # of these candidates, find the most extreme values which mainting the required number of overlapping timepoints to be considered.
  num.overlapping.points <- sapply(neg.extreme.candidates, FUN=calc_num_overlapping_points, original=original)
  if (all(num.overlapping.points < min.num.overlapping.points)) {
    stop(paste0('calc_extreme_shifts():\nafter applying stretch factor:', stretch, ' to ', transformed.timecourse, ', none of the considered shifts have ',
                'min.num.overlapping.points (', min.num.overlapping.points, ') overlapping timepoints with the other timecourse!\n',
                "maybe try a smaller stretch, and double check you're applying it to the correct timecourse." ))
  }

  neg.extreme <- min(neg.extreme.candidates[num.overlapping.points >= min.num.overlapping.points])

  num.overlapping.points <- sapply(pos.extreme.candidates, FUN=calc_num_overlapping_points, original=original)
  pos.extreme <- max(pos.extreme.candidates[num.overlapping.points >= min.num.overlapping.points])

  # hard code maximum and minimum allowed shifts, as noticed spurious registrations when too extreme shifts
  # allowed
  if (neg.extreme < (-1*shift.extreme)) {
    neg.extreme <- -1 * shift.extreme
  }
  if (pos.extreme > 1*shift.extreme) {
    pos.extreme <- shift.extreme
  }

  return(list(neg.extreme, pos.extreme))
}



#' @export
apply_stretch <- function(mean.df, best_shifts) {
  # gets the applied stretch from the best_shifts df
  test <- data.table::copy(mean.df)

  # get the stretch factor applied in the best_shifts
  #stopifnot(length(unique(best_shifts$stretch))==1)
  #stretch_factor <- unique(best_shifts$stretch)

  #stopifnot(length(unique(mean.df$locus_name))==length(best_shifts$gene)) # best_shifts should have 1 row per gene
  # if (length(unique(mean.df$locus_name)) != length(best_shifts$gene)) {
  #   print('ERROR!:')
  #   print(length(unique(mean.df$locus_name)))
  #   print(length(best_shifts$gene))
  #   stop()
  # }


  # stretch the arabidopsis expression data, leave the rapa as is
  test[, delta.time:=timepoint - min(timepoint), by=.(accession)]
  ro18.test <- test[test$accession=='Ro18',]
  col0.test <- test[test$accession=='Col0']
  #test$delta.time[test$accession=='Col0'] <- test$delta.time[test$accession=='Col0']*stretch_factor
  col0.test <- merge(col0.test, best_shifts[, c('gene', 'stretch')], by.x='locus_name', by.y='gene')
  col0.test$delta.time <- col0.test$delta.time * col0.test$stretch
  col0.test$stretch <- NULL
  test <- rbind(ro18.test, col0.test)

  # record the stretched times (before indiv shifting applied)
  test$stretched.time.delta <- test$delta.time # record the time (from start of timecourse) after stretching,
  test$shifted.time <- test$delta.time
  # after stretching, add the time to the first datapoint (7d for ara, 11d for ro18) back on
  test$shifted.time[test$accession=='Col0'] <- test$shifted.time[test$accession=='Col0'] + 11 #7
  test$shifted.time[test$accession=='Ro18'] <- test$shifted.time[test$accession=='Ro18'] + 11
  test$delta.time <- NULL

  return(test)
}


# mean.df <- all.data.df
# best_shifts
#' @export
apply_best_shift <- function(mean.df, best_shifts) {
  # take unregistered expression over time, and the best shifts, and
  # return the registered expression over time for each gene

  test <- data.table::copy(mean.df)

  test <- apply_stretch(mean.df, best_shifts)

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
    test$shifted.time[test$accession=='Col0' & test$locus_name==curr.gene] <- test$shifted.time[test$accession=='Col0' & test$locus_name==curr.gene] + curr.best.shift

    # tmp <- test[test$locus_name==curr.gene]
    # ggplot2::ggplot(tmp)+
    #   ggplot2::aes(x=shifted.time, y=mean.cpm, color=accession) +
    #   ggplot2::geom_point()
  }

  print('done!')

  return(test)
}


#' @export
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

# test <- curr.data.df
#' @export
get_compared_timepoints <- function(test) {
  # flag the arabidopsis timepoints which overlap the brassica timecourse, and so will be compared
  bra.min <- min(test$shifted.time[test$accession=='Ro18'])
  bra.max <- max(test$shifted.time[test$accession=='Ro18'])

  # get the arabidopsis times which used
  test$is.compared <- FALSE
  test$is.compared[(test$accession=='Col0' & (test$shifted.time >= bra.min & test$shifted.time <=bra.max))] <- TRUE

  # get the extreme brassica times which used - bigger or equal than Ara max, and smaller or equal than Ara min, because have to project
  #  Ara onto Bra
  ara.max <- max(test$shifted.time[test$accession=='Col0' & test$is.compared==TRUE])
  ara.min <- min(test$shifted.time[test$accession=='Col0' & test$is.compared==TRUE])
  bra.max <- max_is_compared_to_arabidopsis(ara.max, test[test$accession=='Ro18', ])
  bra.min <- min_is_compared_to_arabidopsis(ara.min, test[test$accession=='Ro18', ])

  # use these to get all the brassica times which used
  test$is.compared[(test$accession=='Ro18' & (test$shifted.time >= bra.min & test$shifted.time <=bra.max))] <- TRUE

  return(test)
}

#ara.expression <- ara.compared$mean.cpm
#bra.expression <- ara.compared$pred.bra.expression
#' @export
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
#' @export
max_is_compared_to_arabidopsis <- function(arabidopsis.time, bra.dt) {
  # return the largest brassica time which is used in comparison t the arabidopsis time
  # the smallest one greater to or equal to arabidopsis time

  # if using for rep data, then repeats of the same points screws it up
  bra.dt <- unique(subset(bra.dt, select=c('timepoint', 'shifted.time')))

  # return the bra dt shifted timepoint which is greater than, or equal to the Ara time.
  bra.dt$diff <- bra.dt$shifted.time - arabidopsis.time
  candidates <- bra.dt[bra.dt$diff>=0,]
  bra.max.time <- candidates$shifted.time[candidates$diff==min(candidates$diff)]

  return(bra.max.time)

  #data.table::setorder(bra.dt, diff)
  #nearest.points <- bra.dt[1:2,]
  #data.table::setorder(nearest.points, shifted.time)
  #return(nearest.points$shifted.time[2])
}

#arabidopsis.time=ara.min
#' @export
min_is_compared_to_arabidopsis <- function(arabidopsis.time, bra.dt) {
  # return the smallest brassica time which is used in comparison t the arabidopsis time
  # the biggest one smaller than or equal to the arabidopsis time

  # if using for rep data, then repeats of the same points screws it up
  bra.dt <- unique(subset(bra.dt, select=c('timepoint', 'shifted.time')))

  bra.dt$diff <- bra.dt$shifted.time - arabidopsis.time
  candidates <- bra.dt[bra.dt$diff<=0,]
  bra.min.time <- candidates$shifted.time[candidates$diff==max(candidates$diff)]
  return(bra.min.time)

}

# arabidopsis.time <- 2
# bra.dt <- ara.df
#' @export
interpolate_brassica_comparison_expression <- function(arabidopsis.time, bra.dt) {

  # arabidopsis time is outside of the range of the bra.dt shifted timepoints
  bra.dt$diff <- bra.dt$shifted.time - arabidopsis.time

  # if outside of comparible range (time is smaller than all bra.dt time or bigger than all)
  if (all(bra.dt$diff > 0) | all(bra.dt$diff < 0) ) {
    return(NA)
  }

  # otherwise,  cut down brassica observations to the two nearest timepoints to the arabidopsis time
  bra.dt$diff <- abs(bra.dt$shifted.time - arabidopsis.time)
  data.table::setorder(bra.dt, diff)
  nearest.points <- bra.dt[1:2,]

  # linearly interpolate between these points to estimate the comparison expression value
  data.table::setorder(nearest.points, shifted.time) # so [1] is earlier time
  time.diff <- nearest.points$shifted.time[2] - nearest.points$shifted.time[1] #
  expression.diff <- nearest.points$mean.cpm[2] - nearest.points$mean.cpm[1]
  grad <- expression.diff / time.diff
  pred.expression <- nearest.points$mean.cpm[1] + (nearest.points$diff[1]) * grad

  return(pred.expression)
}


#get_shifted_expression(shift_results, exp)
#' @export
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
