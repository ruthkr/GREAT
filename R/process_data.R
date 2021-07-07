# stretches=c(1, 1.5, 2.0)
# initial.rescale=FALSE
# do_rescale=FALSE
# min_num_overlapping_points = 4
# test.genes <- unique(mean_df$locus_name)[1:101]
# #test.genes <- 'MSTRG.10244'
# mean_df <- mean_df[mean_df$locus_name %in% test.genes,]
# all_data_df <- all_data_df[all_data_df$locus_name %in% test.genes,]
# shift_extreme=4
# transformed.timecourse <- 'Ro18'

#' @export
prepare_scaled_and_registered_data <- function(mean_df, all_data_df, stretches, initial.rescale, do_rescale, min_num_overlapping_points, shift_extreme, transformed.timecourse) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  ## APPLY NORMALISATION OF EXPRESSION FOR EACH GENE ACROSS ALL TIMEPOINTS ##

  # hardcoded all the functions to use 'Col0', and 'Ro18'. Rather than fix all instances
  # of that, just temporarily rename them here, and turn back at the end.
  # will apply stretch, and shift to the "transformed.timecourse" accession (which should be the quicker one...)
  L <- change.accession.names(mean_df, all_data_df, transformed.timecourse)
  mean_df <- L[['mean_df']]
  all_data_df <- L[['all_data_df']]
  original.transformed.accession <- L[['original.transformed.accession.name']]
  original.other.accession <- L[['original.other.accession.name']]

  mean_df.sc <- data.table::copy(mean_df)
  # specify what kind of scaling
  # TODO: handle my_scale in conditional
  mean_df.sc[, sc.mean_cpm:=scale(mean_cpm, scale=TRUE, center=TRUE), by=.(locus_name, accession)]
  #mean_df.sc[, sc.mean_cpm:=my_scale(mean_cpm), by=.(locus_name, accession)]

  ## APPLY INDIVIDUAL SHIFTING ---
  # optimise transformations applied to arabidopsis gene profile to map onto the brassicas - shift in x direction, using mean for mapping,
  # and only rescale expression using the candidate timepoints in common.

  # specify which data to use for registering.
  # whether prior rescaled mean, or mean data should be used for registration
  if (initial.rescale==TRUE) {
    # apply rescale to mean_df prior to registration
    to_shift_df <- data.table::copy(mean_df.sc)
    to_shift_df$mean_cpm <- to_shift_df$sc.mean_cpm
    to_shift_df$sc.mean_cpm <- NULL

    # apply THE SAME rescale to all_data_df prior to registration
    # TODO: handle my_scale in conditional
    all_data_df <- scale_all_rep_data(mean_df, all_data_df, 'scale')
    #all_data_df <- scale_all_rep_data(mean_df, all_data_df, 'my_scale')


    # sanity plot that rescale all data worked
    # ggplot2::ggplot(all_data_df[all_data_df$locus_name=='BRAA01G000040.3C'])+
    #   ggplot2::aes(x=timepoint, y=mean_cpm, color=accession)
    #   ggplot2::geom_point()
  } else {
    to_shift_df <- data.table::copy(mean_df)
  }

  print(paste0('Max value of mean_cpm of all_data_df :', max(all_data_df$mean_cpm)))
  # ggplot2::ggplot(to_shift_df[to_shift_df$locus_name=='BRAA03G004600.3C'])+
  #   ggplot2::aes(x=timepoint, y=mean_cpm, color=accession)
  #   ggplot2::geom_point()
  # tst <- all_data_df
  # ggplot2::ggplot(tst[tst$locus_name=='BRAA01G000040.3C'])+
  #   ggplot2::aes(x=timepoint, y=mean_cpm, color=accession) +
  #   ggplot2::geom_point()

  # calculate the best registration. Returns all tried registrations, best stretch and shift combo,
  # and AIC/BIC stats for comparison of best registration model to seperate models for expression of
  # each gene in Ro18 and Col0.
  L <- get_best_stretch_and_shift(to_shift_df, all_data_df, stretches, do_rescale, min_num_overlapping_points, shift_extreme)
  all_shifts <- L[['all_shifts']]
  best_shifts <- L[['best_shifts']]
  model_comparison_dt <- L[['model_comparison_dt']]

  print(paste0('Max value of all_shifts mean_cpm :', max(all_shifts$mean_cpm)))


  # report model comparison results
  model_comparison_dt$BIC.registered.is.better <- (model_comparison_dt$registered.BIC < model_comparison_dt$seperate.BIC)
  model_comparison_dt$AIC.registered.is.better <- (model_comparison_dt$registered.AIC < model_comparison_dt$seperate.AIC)
  model_comparison_dt$ABIC.registered.is.better <- (model_comparison_dt$BIC.registered.is.better & model_comparison_dt$AIC.registered.is.better)
  print('################## Model comparison results #######################')
  print(paste0('AIC finds registration better than seperate for :', sum(model_comparison_dt$AIC.registered.is.better), ' / ', nrow(model_comparison_dt)))
  print(paste0('BIC finds registration better than seperate for :', sum(model_comparison_dt$BIC.registered.is.better), ' / ', nrow(model_comparison_dt)))
  print(paste0('AIC & BIC finds registration better than seperate for :', sum(model_comparison_dt$ABIC.registered.is.better), ' / ', nrow(model_comparison_dt)))
  print('###################################################################')


  # get the best-shifted and stretched mean gene expression, only to genes which registration is better than
  # seperate models by BIC. Don't stretch out, or shift genes for which seperate is better.
  # registration is applied to col0.
  shifted.mean_df <- apply_shift_to_registered_genes_only(to_shift_df, best_shifts, model_comparison_dt)
  # shifted.mean_df <- apply_shift_to_all(to_shift_df, best_shifts, model_comparison_dt)
  print(paste0('Max value of mean_cpm :', max(shifted.mean_df$mean_cpm)))
  # shifted.mean_df <- apply_best_shift(to_shift_df, best_shifts) # can be NA if exactly tied for what the best shift was

  # GOI <- 'MSTRG.11237'
  # ggplot2::ggplot(all_data_df[all_data_df$locus_name==GOI])+
  #   ggplot2::aes(x=timepoint, y= mean_cpm, color=accession) +
  #   ggplot2::geom_point()
  #
  # #sanity plot that done right
  # ggplot2::ggplot(shifted.mean_df[shifted.mean_df$locus_name==GOI])+
  #   ggplot2::aes(x=shifted_time, y=mean_cpm, color=accession) +
  #   ggplot2::geom_point()+
  #   ggplot2::geom_line()


  # impute arabidopsis values at times == to the observed brassica points for each shifted arabidopsis gene
  # so can compare using heatmap.
  # arabidopsis curves are the ones that been shifted around. Linear impute values for these
  # curves so that brassica samples can be compared to an arabidopsis point.
  imputed.mean_df <- impute_arabidopsis_values(shifted.mean_df)

  #sanity plot that done right
  # ggplot2::ggplot(shifted.mean_df[shifted.mean_df$locus_name=='BRAA01G001540.3C'])+
  #   ggplot2::aes(x=shifted_time, y=mean_cpm, color=accession) +
  #   ggplot2::geom_point()+
  #   ggplot2::geom_line()
  # ggplot2::ggplot(imputed.mean_df[imputed.mean_df$locus_name=='BRAA01G001540.3C'])+
  #   ggplot2::aes(x=shifted_time, y=mean_cpm, color=accession)+
  #   ggplot2::geom_point()+
  #   ggplot2::geom_line()


  # fix the accession names to the ones actually passed in:
  # mean_df <- fix.accessions(mean_df, original.transformed.accession, original.other.accession)
  # mean_df.sc <- fix.accessions(mean_df.sc, original.transformed.accession, original.other.accession)
  # imputed.mean_df <- fix.accessions(imputed.mean_df, original.transformed.accession, original.other.accession)


  OUT <- list('mean_df'=mean_df,
              'mean_df.sc'=mean_df.sc,
              'imputed.mean_df'=imputed.mean_df,
              'all_shifts'=all_shifts,
              'model.comparison'=model_comparison_dt)
}



#' @export
change.accession.names <- function(mean_df, all_data_df, transformed.timecourse) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
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

  # change all_data_df
  new.all_data_df.accession <- all_data_df$accession
  new.all_data_df.accession[all_data_df$accession==transformed.timecourse] <- 'Col0'
  new.all_data_df.accession[all_data_df$accession!=transformed.timecourse] <- 'Ro18'
  all_data_df$accession <- new.all_data_df.accession

  return(list('mean_df'=mean_df,
              'all_data_df'=all_data_df,
              'original.transformed.accession.name'=original.transformed.timecourse.name,
              'original.other.accession.name'=original.other.accession.name))
}

#all.rep.data <- all_data_df
#' @export
scale_all_rep_data <- function(mean_df, all.rep.data, scale.func) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # apply the same scaling which done to the mean expression data
  # to all the reps.
  # (subtract mean, and divide by sd), using the values for the mean data
  # as this is what was used to find the best shift.

  # calculate the summary stats to use for the rescaling
  gene.expression.stats <- unique(mean_df[,
                                          .(mean_val=mean(mean_cpm),
                                            sd_val=stats::sd(mean_cpm)),
                                          by=.(locus_name, accession)])

  all.rep.data <- merge(all.rep.data, gene.expression.stats, by=c('locus_name', 'accession'))
  if (scale.func == 'scale') {
    all.rep.data$scaled.norm.cpm <- (all.rep.data$mean_cpm - all.rep.data$mean_val) / all.rep.data$sd_val
  } else if (scale.func == 'my_scale') {
    all.rep.data$scaled.norm.cpm <- (all.rep.data$mean_cpm / all.rep.data$mean_val)
  } else {
    print('invalid scale option for scale_all_rep_data')
    stop()
  }

  out <- subset(all.rep.data, select=c('locus_name', 'accession', 'tissue', 'timepoint',
                                       'scaled.norm.cpm'))

  names(out)[names(out)=='scaled.norm.cpm'] <- 'mean_cpm'

  # ggplot2::ggplot(mean_df[mean_df$locus_name=='BRAA01G000040.3C', ], ) +
  #  ggplot2::aes(x=timepoint, y=mean_cpm, color=accession
  #  ggplot2::geom_point()


  return(out)
}


#' Calculate best shifts and stretches for each gene, also calculate AIC/BIC under registration or non-registration
#'
#' `get_best_stretch_and_shift` is a function to stretch in all stretches and calculates best shift, by comparing SUM of squares difference. For the best shift in each stretch, compares to separate models to calculate AIC/BIC under registration or no registration.
#'
#' @param to_shift_df Input data containing mean of each time point.
#' @param all_data_df Input all data (without taking mean).
#' @param stretches Vector data of stretches.
#' @param do_rescale Apply "scale" to compared points for each shift if TRUE, use original mean expression data if FALSE.
#' @param min_num_overlapping_points Bound the extreme allowed shifts, such than at least this many timepoints are being compared for both accessions.
#' @param shift_extreme Approximation of maximum and minimum shifts allowed.
#' @param num_shifts Number of different shifts to be considered.
#' @param testing Showing a plot of the progress if TRUE, otherwise if FALSE.
#' @param accession_data_to_align Accession name of data which will be aligned.
#' @param accession_data_target Accession name of data target.
#' @param data_to_align_time_added Time points to be added in data to align.
#' @param data_target_time_added Time points to be added in data target.
#'
#' @return List of data frames (a) all_shifts : all the combos of stretching and shifting tried for each gene, (b) best_shifts : the best stretch and shift combo found for each gene, as well as info for scaling, and (c) model_comparison.dt : AIC / BIC scores for best registerd model found, compared to seperate model for each genes expression in the 2 accessions.
#' @export
get_best_stretch_and_shift <- function(to_shift_df,
                                       all_data_df,
                                       stretches,
                                       do_rescale,
                                       min_num_overlapping_points,
                                       shift_extreme,
                                       num_shifts,
                                       testing,
                                       accession_data_to_align,
                                       accession_data_target,
                                       data_to_align_time_added,
                                       data_target_time_added) {

  # Warning to make sure users have correct accession data
  if (!('Col0' %in% all_data_df$accession & 'Ro18' %in% all_data_df$accession)) {
    stop('get_best_stretch_and_shift() : data accessions should have been
         converted to correct accession.')
  }

  all_all_shifts <- rep(list(0), length(stretches))
  all_best_shifts <- rep(list(0), length(stretches))
  all_model_comparison_dt <- rep(list(0), length(stretches))


  for (i in 1:length(stretches)) {
    stretch <- stretches[i]
    message(paste0('testing models for stretch factor = ', stretch))

    # Calculate all the shift scores given this stretch. Score is mean(dist^2), over overlapping points
    # if do_rescale=T, is rescaled by the mean FOR THE OVERLAPPING POINTS. (but not by the SD.)
    all_shifts <- calculate_all_best_shifts(num_shifts,
      mean_df = to_shift_df,
      stretch_factor = stretch,
      do_rescale,
      shift_extreme,
      min_num_overlapping_points,
      testing = FALSE,
      accession_data_to_align,
      accession_data_target)

    all_shifts <- unique(all_shifts) # ensure no duplicated rows

    # Cut down to single best shift for each gene
    all_shifts[, is_best:=get_best_result(.SD), by=.(gene)]
    best_shifts <- all_shifts[is_best==TRUE,]
    all_shifts$is_best <- NULL

    if (nrow(best_shifts)!= length(unique(all_data_df$locus_name))) {
      stop('get_best_stretch_and_shift() : got non-unique best shifts in best_shifts')
    }

    # Calculate the BIC & AIC for the best shifts found with this stretch.compared to treating the
    # gene's expression separately in data to align and data target
    model_comparison_dt <- calculate_all_model_comparison_stats(all_data_df,
                                                                best_shifts,
                                                                accession_data_to_align,
                                                                accession_data_target,
                                                                data_to_align_time_added,
                                                                data_target_time_added)


    # Add info on the stretch and shift applied
    model_comparison_dt <- merge(model_comparison_dt, best_shifts[, c('gene', 'stretch', 'shift'),],
                                 by='gene')

    # Record the results for the current stretch factor
    all_all_shifts[[i]] <- all_shifts
    all_best_shifts[[i]] <- best_shifts
    all_model_comparison_dt[[i]] <- model_comparison_dt
    message(paste0('finished testing models for stretch factor = ', stretch), "\n")

  }

  all_shifts <- do.call('rbind', all_all_shifts) # all the combinations of shift, and stretch tried
  all_best_shifts <- do.call('rbind', all_best_shifts) # the best shifts for each stretch
  all_model_comparison_dt <- do.call('rbind', all_model_comparison_dt) # model comparison of best shift (for each stretch) to seperate modeles

  # Get the best registration applied (best stretch, and best shift) for each gene,
  # picking by BIC alone will favour fewer overlapping (considered) data points.
  # Pick best in order to maximise how much better register.BIC is than separate.BIC
  all_model_comparison_dt$delta.BIC <- all_model_comparison_dt$registered.BIC - all_model_comparison_dt$seperate.BIC

  # Best is one for which registered.BIC is as small as possible compared to separate.BIC
  all_model_comparison_dt[, is_best := (delta.BIC == min(delta.BIC)), by = .(gene)]
  best_model_comparison.dt <- all_model_comparison_dt[all_model_comparison_dt$is_best==TRUE]

  # If there is a tie for best registration for a gene, keep the first one as the best
  if (any(duplicated(best_model_comparison.dt$gene))) {
    print(paste0('found ', sum(duplicated(best_model_comparison.dt$gene)), ' tied optimal registrations. Removing dupliates'))
    best_model_comparison.dt <- best_model_comparison.dt[!(duplicated(best_model_comparison.dt$gene)),]
  }

  best_model_comparison.dt$delta.BIC <- NULL

  # Cut down best shifts to the best shift for the best stretch only
  best_shifts <- merge(all_best_shifts,
                       best_model_comparison.dt[, c('gene', 'stretch', 'shift')],
                       by=c('gene', 'stretch', 'shift'))

  # There should be only 1 best shift for each gene, stop if it is not the case
  stopifnot(nrow(best_shifts) == length(unique(to_shift_df$locus_name)))

  return(list('all_shifts'=all_shifts,
              'best_shifts'=best_shifts,
              'model_comparison_dt'=best_model_comparison.dt))

}

#' @export
apply_shift_to_registered_genes_only <- function(to_shift_df,
                                                 best_shifts,
                                                 model_comparison_dt) {

  # genes for which registration model is better than seperate model
  genes.to.register <- model_comparison_dt$gene[model_comparison_dt$BIC.registered.is.better]
  # apply the registration transformation to these genes
  if (length(genes.to.register > 0)) {
    register.dt <- to_shift_df[to_shift_df$locus_name %in% genes.to.register,]
    registered.dt <- apply_best_shift(register.dt, best_shifts)
    registered.dt$is.registered <- TRUE
  }

  # genes for which the seperate model is better than registration model
  genes.to.keep.seperate <- model_comparison_dt$gene[!(model_comparison_dt$BIC.registered.is.better)]

  # generate the columns for these needed to concat. with registered.dt
  seperate.dt <- to_shift_df[to_shift_df$locus_name %in% genes.to.keep.seperate,]
  seperate.dt$stretched.time.delta <- 0 # in order to ensure that seperate copy
  # print('line 594')
  # print(min(timepoint))
  seperate.dt[, stretched.time.delta:=timepoint - min(timepoint), by=.(locus_name, accession)]
  seperate.dt$shifted_time <- seperate.dt$stretched.time.delta + 14 # add eleven, as this is done for the registered genes
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
apply_shift_to_all <- function(to_shift_df, best_shifts, model_comparison_dt) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])

  registered.dt <- to_shift_df
  registered.dt <- apply_best_shift(registered.dt, best_shifts)
  # registered.dt$shifted_time <- registered.dt$stretched.time.delta + 14 # add eleven, as this is done for the registered genes

  registered.dt$is.registered <- TRUE
  return(registered.dt)
}

#' @export
impute_arabidopsis_values <- function(shifted.mean_df) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # Arabidopsis gene expression profiles are shifted all over. Need to impute times at set of common timepoints
  # in order to allow sample distance comparison to Ro18.
  #
  # Ro18 genes haven't been shifted around, therefore imputed timepoints are realtive to Ro18 timepoints.
  # We ONLY have Ro18 observations for 11, 13, 15, ... 35
  #
  # Col0 can have had different shifts applied to them, so need to impute them to a common scale, and to compare to Ro18
  #
  # therefore, we only need to impute Arabidopsis gene expression, which will be compared to the Ro18 timepoints.
  # BUT don't want to discard any col0 data, so want to generate imputed time observations for col0 from minimum
  # to maximum shifted timepoints for Col0 not just for 11,13,15, etc Ro18 observations.


  # sanity plotting - a registered one
  # ggplot2::ggplot(shifted.mean_df[shifted.mean_df$locus_name=='BRAA01G001090.3C', ])+
  #   ggplot2::aes(x=shifted_time, y=mean_cpm, color=accession)+
  #   ggplot2::geom_point()
  # # an unregistered one
  # ggplot2::ggplot(shifted.mean_df[shifted.mean_df$locus_name=='BRAA01G003140.3C', ])+
  #   ggplot2::aes(x=shifted_time, y=mean_cpm, color=accession)+
  #   ggplot2::geom_point()

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
    # if (all(is.na(curr.df$mean_cpm))) {
    #   next
    # }

    ara.df <- curr.df[curr.df$accession=='Col0',]
    bra.df <- curr.df[curr.df$accession=='Ro18',]

    interp.ara.df <- data.table::data.table(data.frame('locus_name'=curr.gene, 'accession'='Col0', 'tissue'='apex', 'timepoint'=NA,
                                                       'stretched.time.delta'= NA, 'shifted_time'=imputed.timepoints,
                                                       'is.registered'= unique(ara.df$is.registered)[1]))

    # for each brassica timepoint, interpolate the comparible arabidopsis expression
    # by linear interpolation between the neighbouring 2 ara values. If not between 2 ara values
    # because shifted outside comparible range, set to NA
    interp.ara.df$mean_cpm <- sapply(imputed.timepoints, interpolate_data_target_comparison_expression, bra.dt=ara.df)


    # sanity testing - line is interpolated.
    # ggplot2::ggplot(curr.df)+
    #   ggplot2::aes(x=shifted_time, y=mean_cpm, color=accession)+
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
