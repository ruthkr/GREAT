
#' @export
get_best_result <- function(df) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
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

#' @export
calculate_all_model_comparison_stats <- function(all.data.df, best_shifts) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
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
  #   ggplot2::aes(x=shifted_time, y=mean.cpm, color=accession) +
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

# mean_df <- all.data.df
# best_shifts
#' @export
apply_best_shift <- function(mean_df, best_shifts) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # take unregistered expression over time, and the best shifts, and
  # return the registered expression over time for each gene

  test <- data.table::copy(mean_df)

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
  # curr.gene <- 'BRAA01G000040.3C'
  #curr.gene <- unique(test$locus_name)[1]
  for (curr.gene in unique(test$locus_name)) {
    #print(curr.gene)
    curr.best.shift <- best_shifts$shift[best_shifts$gene==curr.gene]
    test$shifted_time[test$accession=='Col0' & test$locus_name==curr.gene] <- test$shifted_time[test$accession=='Col0' & test$locus_name==curr.gene] + curr.best.shift

    # tmp <- test[test$locus_name==curr.gene]
    # ggplot2::ggplot(tmp)+
    #   ggplot2::aes(x=shifted_time, y=mean.cpm, color=accession) +
    #   ggplot2::geom_point()
  }

  print('done!')

  return(test)
}

#' @export
apply_stretch <- function(mean_df, best_shifts) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # gets the applied stretch from the best_shifts df
  test <- data.table::copy(mean_df)

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
  Ro18.test <- test[test$accession=='Ro18',]
  col0.test <- test[test$accession=='Col0']
  #test$delta_time[test$accession=='Col0'] <- test$delta_time[test$accession=='Col0']*stretch_factor
  col0.test <- merge(col0.test, best_shifts[, c('gene', 'stretch')], by.x='locus_name', by.y='gene')
  col0.test$delta_time <- col0.test$delta_time * col0.test$stretch
  col0.test$stretch <- NULL
  test <- rbind(Ro18.test, col0.test)

  # record the stretched times (before indiv shifting applied)
  test$stretched.time.delta <- test$delta_time # record the time (from start of timecourse) after stretching,
  test$shifted_time <- test$delta_time
  # after stretching, add the time to the first datapoint (7d for ara, 11d for Ro18) back on
  test$shifted_time[test$accession=='Col0'] <- test$shifted_time[test$accession=='Col0'] + 14
  test$shifted_time[test$accession=='Ro18'] <- test$shifted_time[test$accession=='Ro18'] + 14
  test$delta_time <- NULL

  return(test)
}

#' @export
apply_best_normalisation <- function(test, best_shifts) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
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

# curr.sym <- 'BRAA01G001320.3C'
# all.data.df <- all.data.df
# is.testing <- TRUE
#' @export
compare_registered_to_unregistered_model <- function(curr.sym, all.data.df, is.testing) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # compare the overlapping timepoints in brassica and arabidopsis after the best registration,
  # and without registration (use the same timepoints for both models).
  # use the stretched data for both models, whether considering as registered or not.


  curr.data.df <- all.data.df[all.data.df$locus_name==curr.sym]

  # print('line 662')
  # print(curr.data.df)

  # flag the timepoints to be used in the modelling, only the ones which overlap!
  curr.data.df <- get_compared_timepoints(curr.data.df)

  # ggplot2::ggplot(curr.data.df)+
  #   ggplot2::aes(x=shifted_time, y=mean.cpm, shape=is.compared, color=accession)+
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


  ara.fit <- stats::lm(mean.cpm~splines::bs(shifted_time, df=num.spline.params, degree=3), data=ara.spline.data)
  bra.fit <- stats::lm(mean.cpm~splines::bs(shifted_time, df=num.spline.params, degree=3), data=bra.spline.data)
  combined.fit <- stats::lm(mean.cpm~splines::bs(shifted_time, df=num.spline.params, degree=3), data=combined.spline.data)
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
    ara.pred.df <- unique(data.frame('shifted_time'=ara.spline.data$shifted_time,
                                     'mean.cpm'=ara.pred, 'accession'='Col0'))
    bra.pred <- stats::predict(bra.fit)
    bra.pred.df <- unique(data.frame('shifted_time'=bra.spline.data$shifted_time,
                                     'mean.cpm'=bra.pred, 'accession'='Ro18'))

    combined.pred <- stats::predict(combined.fit)
    combined.pred.df <- unique(data.frame('shifted_time'=combined.spline.data$shifted_time,
                                          'mean.cpm'=combined.pred, 'accession'='registered'))
    spline.df <- rbind(ara.pred.df, bra.pred.df, combined.pred.df)

    p <- ggplot2::ggplot(data=combined.spline.data)+
      ggplot2::aes(x=shifted_time, y=mean.cpm, colour=accession) +
      ggplot2::geom_point()+
      ggplot2::geom_line(data=spline.df)+
      ggplot2::ggtitle(paste0(curr.sym, ' : sep AIC:combo AIC=', round(seperate.AIC), ':', round(combined.AIC),
                              ', sep BIC: combo BIC=', round(seperate.BIC), ':', round(combined.BIC)))
    # ggplot2::ggsave(paste0('./testing/fitted_splines/', curr.sym, '_', max(ara.pred.df$shifted_time), '.pdf'))
  }

  return(list(seperate.AIC, combined.AIC, seperate.BIC,combined.BIC))
}



# for all genes, in mean_df, get the scores, and the shifts which result in them after
# stretching arabidopsis by "stretch_factor, and applying shifts forward and
# backward, with the extremes defined by allowing 5 overlapping points for comparison.
# mean_df <- to.shift.df
# stretch_factor <- 2
# do_rescale=F
#' @export
calculate_all_best_shifts <- function(mean_df, stretch_factor, do_rescale, min_num_overlapping_points, shift_extreme) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])

  # Initialize vectors
  symbols <- c()
  num_points <- c()
  all_scores_list <- rep(list(0), length(unique(mean_df$locus_name)))

  # Get the extreme shifts which can be applied to the genes
  M <- get_extreme_shifts_for_all(mean_df, stretch_factor, min_num_overlapping_points, shift_extreme)
  min_shift <- M[[1]]
  max_shift <- M[[2]]
  print(paste0("min shift:", min_shift, "max shift:", max_shift))


  count <- 0
  for (i in 1:length(unique(mean_df$locus_name))) {
    #for (curr_sym in unique(mean_df$locus_name)) {
    curr_sym <- unique(mean_df$locus_name)[i]
    if (count %% 100 == 0) {
      print(paste0(count, ' / ', length(unique(mean_df$locus_name))))
    }

    # out is mean SSD between arabidopsis, and interpolated brassica (interpolated between 2 nearest points)
    # ggplot2::ggplot(mean_df[mean_df$locus_name==curr_sym,])+
    #   ggplot2::aes(x=timepoint, y=mean.cpm, color=accession) +
    #   ggplot2::geom_point()

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

    count <- count + 1
  }

  all_scores_df <- do.call('rbind', all_scores_list)

  return(all_scores_df)
}

# test <- mean_df
# stretch_factor
# min_num_overlapping_points
# shift_extreme
#' @export
get_extreme_shifts_for_all <- function(test, stretch_factor, min_num_overlapping_points, shift_extreme) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
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
  #data.table::setorder(test, delta_time)
  #min_shift <- min(test$delta_time[test$accession=='Ro18']) - test$delta_time[test$accession=='Col0'][length(test$delta_time[test$accession=='Col0'])-4]
  #max_shift <- max(test$delta_time[test$accession=='Ro18']) - test$delta_time[test$accession=='Col0'][5]

  M <- calc_extreme_shifts(test, min_num_overlapping_points, shift_extreme)
  return(M)
}

#' @export
calc_extreme_shifts <- function(test, min_num_overlapping_points, shift_extreme) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # calculate the minimum and maximum shifts can apply to Col-0 after the stretch transformation, whilst
  # preserving the criteria that at least min_num_overlapping_points are being compared from both accessions.

  original <- data.table::copy(test)
  original$shifted_time <- original$delta_time

  # print('line 1803')
  # print(original)

  # -ve extreme shift will be -1*exactly the difference between 1 of the stretched Col0 timepoints, and the smallest Ro18 timepoint
  # +ve extreme will be the difference between 1 of the col0 timepoints, and the maximum Ro18 timepoint
  neg_extreme_candidate <- -1*(original$delta_time[original$accession=='Col0'] - min(original$delta_time[original$accession=='Ro18']))
  pos.extreme.candidates <- max(original$delta_time[original$accession=='Ro18']) - original$delta_time[original$accession=='Col0']

  # of these candidates, find the most extreme values which mainting the required number of overlapping timepoints to be considered.
  num.overlapping.points <- sapply(neg_extreme_candidate, FUN=calc_num_overlapping_points, original=original)
  if (all(num.overlapping.points < min_num_overlapping_points)) {
    stop(paste0('calc_extreme_shifts():\nafter applying stretch factor:', stretch, ' to ', transformed.timecourse, ', none of the considered shifts have ',
                'min_num_overlapping_points (', min_num_overlapping_points, ') overlapping timepoints with the other timecourse!\n',
                "maybe try a smaller stretch, and double check you're applying it to the correct timecourse." ))
  }

  neg.extreme <- min(neg_extreme_candidate[num.overlapping.points >= min_num_overlapping_points])

  num.overlapping.points <- sapply(pos.extreme.candidates, FUN=calc_num_overlapping_points, original=original)
  pos.extreme <- max(pos.extreme.candidates[num.overlapping.points >= min_num_overlapping_points])

  # hard code maximum and minimum allowed shifts, as noticed spurious registrations when too extreme shifts
  # allowed
  if (neg.extreme < (-1*shift_extreme)) {
    neg.extreme <- -1 * shift_extreme
  }
  if (pos.extreme > 1*shift_extreme) {
    pos.extreme <- shift_extreme
  }

  return(list(neg.extreme, pos.extreme))
}


# wroking on new implementation of the funciton
# curr_sym <- 'BRAA02G015410.3C'
#test <- mean_df
# stretch_factor <- 1
# do_rescale <- TRUE
#testing=T

#' @export
get_best_shift_new <- function(curr_sym, test, stretch_factor, do_rescale, min_shift, max_shift, testing=FALSE) {
  message_function_header(unlist(stringr::str_split(deparse(sys.call()), "\\("))[[1]])
  # for the current gene, and current stretch_factor, calculate the score for all
  # shifts, and return the scores for all as a table, and the value of the optimal shift.

  # Shift extremes are defined s.t. at least 5 points are compared.

  # do_rescale == TRUE, means apply "scale" to compared points for each shift. ==FALSE, means use original mean expression data

  num.shifts <- 25 # the number of different shifts to be considered.


  test <- test[test$locus_name==curr_sym, ]

  # transform timepoint to be time from first timepoint
  test[, delta_time:=timepoint - min(timepoint), by=.(accession)]
  # apply stretch_factor to the arabidopsis, leave the rapa as is
  test$delta_time[test$accession == 'Col0'] <- test$delta_time[test$accession=='Col0'] * stretch_factor


  # # Try to make it consistent as in apply stretch
  # test$delta_time[test$accession=='Col0'] <- test$delta_time[test$accession=='Col0'] + 7
  # test$delta_time[test$accession=='Ro18'] <- test$delta_time[test$accession=='Ro18'] + 14

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

    # shift the arabidopsis expression timings
    test$shifted_time <- test$delta_time
    test$shifted_time[test$accession == 'Col0'] <- test$delta_time[test$accession == 'Col0'] + curr.shift

    #### test plot - of shifted, UNNORMALISED gene expression
    # if (testing==TRUE) {
    #   p <- ggplot2::ggplot(test)+
    #     ggplot2::aes(x=shifted_time, y=mean.cpm, color=accession)+
    #     ggplot2::geom_point()+
    #     ggplot2::ggtitle(paste0('shift : ', curr.shift))
    #   p
    #   ggplot2::ggsave(paste0('./testing/', curr.shift, '.pdf'))
    # }


    # Cut down to just the arabidopsis and brassica timepoints which compared
    test <- get_compared_timepoints(test)
    compared <- test[test$is.compared==TRUE, ]

    # Renormalise expression using just these timepoints?
    if (do_rescale == TRUE) {
      # record the mean and sd of the compared points, used for rescaling
      # in "apply shift" function
      ara.mean <- mean(compared$mean.cpm[compared$accession=='Col0'])
      bra.mean <- mean(compared$mean.cpm[compared$accession=='Ro18'])
      ara.sd <- stats::sd(compared$mean.cpm[compared$accession=='Col0'])
      bra.sd<- stats::sd(compared$mean.cpm[compared$accession=='Ro18'])

      # do the transformation for here
      if ((ara.sd != 0 | !is.nan(ara.sd)) & (bra.sd != 0 | !is.nan(bra.sd))) { # if neither are 0, so won't be dividing by 0 (which gives NaNs)
        compared[, mean.cpm:=scale(mean.cpm, scale=TRUE, center=TRUE), by=.(accession)]
      } else { # if at least one of them is all 0
        ara.compared <- compared[compared$accession=='Col0',]
        bra.compared <- compared[compared$accession=='Ro18',]
        if((ara.sd == 0) & (bra.sd != 0 | !is.nan(bra.sd))) { # if only ara.sd==0
          bra.compared[, mean.cpm:=scale(mean.cpm, scale=TRUE, center=TRUE), by=.(accession)]
        }
        if ((ara.sd != 0 | !is.nan(ara.sd))  & (bra.sd == 0)) { # if only bra.sd == 0
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
        ggplot2::aes(x=shifted_time, y=mean.cpm, color=accession)+
        ggplot2::geom_point()+
        ggplot2::ggtitle(paste0('shift : ', curr.shift))
      ggplot2::ggsave(paste0('./testing/',stretch_factor, '-', curr.shift, '.pdf'))
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
      print(paste("ara.compared$mean.cpm:", ara.compared$mean.cpm))
      print(ara.compared$pred.bra.expression)
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

  out <- data.table::data.table(data.frame('gene'=curr_sym, 'stretch'=stretch_factor, 'shift'=all.shifts, 'score'=all_scores,
                                           'ara.compared.mean'=all.ara.mean, 'bra.compared.mean'=all.bra.mean,
                                           'ara.compared.sd'=all.ara.sd, 'bra.compared.sd'=all.bra.sd))
  return(out)
}
