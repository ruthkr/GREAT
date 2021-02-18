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
prepare_scaled_and_registered_data <- function(mean.df, all.data.df, stretches, initial.rescale, do.rescale, min.num.overlapping.points, shift.extreme, transformed.timecourse) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
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
  # TODO: handle my.scale in conditional
  mean.df.sc[, sc.mean.cpm:=scale(mean.cpm, scale=TRUE, center=TRUE), by=.(locus_name, accession)]
  #mean.df.sc[, sc.mean.cpm:=my.scale(mean.cpm), by=.(locus_name, accession)]

  ## APPLY INDIVIDUAL SHIFTING ---
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
    # TODO: handle my.scale in conditional
    all.data.df <- scale_all_rep_data(mean.df, all.data.df, 'scale')
    #all.data.df <- scale_all_rep_data(mean.df, all.data.df, 'my.scale')


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
change.accession.names <- function(mean.df, all.data.df, transformed.timecourse) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
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

#all.rep.data <- all.data.df
#' @export
scale_all_rep_data <- function(mean.df, all.rep.data, scale.func) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
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

# to.shift.df
# all.data.df
# stretches
# do.rescale
# min.num.overlapping.points
# shift.extreme
#' @export
get_best_stretch_and_shift <- function(to.shift.df, all.data.df, stretches, do.rescale, min.num.overlapping.points, shift.extreme) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])

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

  # i <- 3 # useless
  for (i in 1:length(stretches)) {
    stretch <- stretches[i]
    message(paste0('testing models for stretch factor = ', stretch))
    # calculate all the shift scores given this stretch. Score is mean(dist^2), over overlapping points
    # if do.rescale=T, is rescaled by the mean FOR THE OVERLAPPING POINTS. (but not by the SD.)

    # ggplot2::ggplot(to.shift.df[to.shift.df$locus_name=='MSTRG.12467',])+
    #   ggplot2::aes(x=timepoint, y=mean.cpm, color=accession) +
    #   ggplot2::geom_point()

    all_shifts <- calculate_all_best_shifts(to.shift.df, stretch_factor=stretch, do.rescale, min.num.overlapping.points, shift.extreme)

    all_shifts <- unique(all_shifts) # ensure no duplicated rows

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
    message(paste0('finished testing models for stretch factor = ', stretch), "\n")
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
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])

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
impute_arabidopsis_values <- function(shifted.mean.df) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
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

#' @export
fix.accessions <- function(df, original.transformed.accession, original.other.accession) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # swap Col0 with original.transformed.accession, and Ro18 with original.other.accession
  new.df.accession <- df$accession
  new.df.accession[df$accession=='Col0'] <- original.transformed.accession
  new.df.accession[df$accession=='Ro18'] <- original.other.accession
  df$accession <- new.df.accession

  return(df)
}
