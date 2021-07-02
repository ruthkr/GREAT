#' Get the best result
#'
#' `get_best_result` is a function to get best result obtained from score calculated from applied shifts and stretch factors.
#'
#' @param df Input data frame containing value of applied shifts, stretches, and calculated score.
#'
#' @return Original data frame input with additional column indicating whether a pair of stretch and shift gives the best score.
#' @export
get_best_result <- function(df) {

  # return TRUE/FALSE vector. TRUE for the smallest score
  # if tied for this, true for the one with the smallest stretch. (1x is smaller than 0.75x though)
  # if tied, then the one with the smallest shift

  is_best <- df$score == min(df$score)

  if(sum(is_best) == 1) {
    return(is_best)
  } else {
    cand_stretches <- df$stretch[is_best]
    # get the stretch with the best score, with the smallest divergence from 1
    min_stretch <- unique(cand_stretches[abs(cand_stretches - 1) == min(abs(cand_stretches - 1))])
    is_best[df$stretch != min_stretch] <- FALSE

    if(sum(is_best) == 1) {
      return(is_best)
    } else {
      cand_shifts <- df$shift[is_best]
      min_shift <- unique(cand_shifts[cand_shifts == min(cand_shifts)])
      is_best[df$shift != min_shift] <- FALSE
      if (sum(is_best) == 1) {
        return(is_best)
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
  #   ggplot2::aes(x=shifted_time, y=mean_cpm, color=accession) +
  #   ggplot2::geom_point()
  # ggplot2::ggplot(tst[tst$locus_name=='BRAA01G000040.3C'])+
  #   ggplot2::aes(x=timepoint, y=mean_cpm, color=accession) +
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
    #   ggplot2::aes(x=shifted_time, y=mean_cpm, color=accession) +
    #   ggplot2::geom_point()
  }

  print('done!')

  return(test)
}


#' Apply stretch factor
#'
#' @param mean_df Input data.
#' @param best_shifts Input dataframe containing information of best shifts.
#' @param accession_data_to_align Accession name of data which will be aligned.
#' @param accession_data_target Accession name of data target.
#' @param data_to_align_time_added Time points to be added in data to align.
#' @param data_target_time_added Time points to be added in data target.
#'
#' @return
#' @export
apply_stretch <- function(mean_df,
                          best_shifts,
                          accession_data_to_align,
                          accession_data_target,
                          data_to_align_time_added = 11,
                          data_target_time_added = 11) {

  data <- data.table::copy(mean_df)

  # Stretch the expression of data to align, leave data target as is
  data[, delta_time := timepoint - min(timepoint), by = .(accession)]

  # Filter data based on the accession
  data_target <- data[data$accession == accession_data_target, ]
  data_to_align <- data[data$accession == accession_data_to_align, ]

  # Get the info of the strecth factor and merge data into one single data frame
  data_to_align <- merge(data_to_align,
    best_shifts[, c("gene", "stretch")],
    by.x = "locus_name",
    by.y = "gene"
  )

  data_to_align$delta_time <- data_to_align$delta_time * data_to_align$stretch
  data_to_align$stretch <- NULL

  # Bind by rows data target and data to align which have been stretched
  data <- rbind(data_target, data_to_align)

  # Record the stretched times (before individual shifting applied)
  data$stretched_time_delta <- data$delta_time # record the time (from start of timecourse) after stretching,
  data$shifted_time <- data$delta_time

  # After stretching, add the time to the first datapoint back on
  data$shifted_time[data$accession == accession_data_to_align] <- data$shifted_time[data$accession == accession_data_to_align] + data_to_align_time_added
  data$shifted_time[data$accession == accession_data_target] <- data$shifted_time[data$accession == accession_data_target] + data_target_time_added
  data$delta_time <- NULL

  return(data)
}


#' @export
apply_best_normalisation <- function(data,
                                     best_shifts,
                                     accession_data_to_align = "Col0",
                                     accession_data_target = "Ro18") {

  # for each gene, in each accession (Ro18 and COl0) normalise by the mean and standard deviation of the compared points.
  # if the gene wasn't compared, set the expresion value to NA

  count = 0
  for (curr.gene in unique(data$locus_name)) {
    if (count %% 100 == 0) {
      print(paste0(count, ' / ', length(unique(data$locus_name))))
    }
    ara.mean <- best_shifts$ara.compared.mean[best_shifts$gene==curr.gene]
    bra.mean <- best_shifts$bra.compared.mean[best_shifts$gene==curr.gene]
    ara.sd<- best_shifts$ara.compared.sd[best_shifts$gene==curr.gene]
    bra.sd<- best_shifts$bra.compared.sd[best_shifts$gene==curr.gene]

    # if was compared
    if (length(ara.mean) != 0) {
      if (ara.sd != 0) { # don't want to divide by 0
        data$mean_cpm[data$locus_name == curr.gene & data$accession == accession_data_to_align] <- (data$mean_cpm[data$locus_name==curr.gene & data$accession == accession_data_to_align] - ara.mean) / ara.sd
      } else {
        data$mean_cpm[data$locus_name == curr.gene & data$accession == accession_data_to_align] <- (data$mean_cpm[data$locus_name == curr.gene & data$accession == accession_data_to_align] - ara.mean)
      }

      if (bra.sd !=0) { # don't want to divide by 0
        data$mean_cpm[data$locus_name==curr.gene & data$accession == accession_data_target] <- (data$mean_cpm[data$locus_name==curr.gene & data$accession == accession_data_target] - bra.mean) / bra.sd
      } else {
        data$mean_cpm[data$locus_name==curr.gene & data$accession == accession_data_target] <- (data$mean_cpm[data$locus_name==curr.gene & data$accession == accession_data_target] - bra.mean)
      }
      if (any(is.na(data$mean_cpm))) {
        print('have NAs in mean_cpm after rescaling in apply best_normalisation() for gene :')
        print(unique(data$locus_name))
        stop()
      }
    } else {
      data$mean_cpm[data$locus_name==curr.gene & data$accession == accession_data_to_align] <- NA
      data$mean_cpm[data$locus_name==curr.gene & data$accession == accession_data_target] <- NA
    }
    count <- count + 1
  }

  return(data)
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
  #   ggplot2::aes(x=shifted_time, y=mean_cpm, shape=is.compared, color=accession)+
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


  ara.fit <- stats::lm(mean_cpm~splines::bs(shifted_time, df=num.spline.params, degree=3), data=ara.spline.data)
  bra.fit <- stats::lm(mean_cpm~splines::bs(shifted_time, df=num.spline.params, degree=3), data=bra.spline.data)
  combined.fit <- stats::lm(mean_cpm~splines::bs(shifted_time, df=num.spline.params, degree=3), data=combined.spline.data)
  # calculate the log likelihoods
  ara.logLik <- stats::logLik(ara.fit)
  bra.logLik <- stats::logLik(bra.fit)
  seperate.logLik <- ara.logLik + bra.logLik # logLikelihoods, so sum
  combined.logLik <- stats::logLik(combined.fit)

  # calculate the comparison.stats - - AIC, BIC, smaller is better!
  # 2*num.spline.params as fitting seperate models for Ara * Col
  seperate.AIC <- calc_AIC(seperate.logLik, 2*num.spline.params)
  combined.AIC <- calc_AIC(combined.logLik, num.spline.params+num.registration.params)

  seperate.BIC <- calc_BIC(seperate.logLik, 2*num.spline.params, num.obs)
  combined.BIC <- calc_BIC(combined.logLik, num.spline.params+num.registration.params, num.obs)


  if (is.testing==TRUE) {
    ara.pred <- stats::predict(ara.fit)
    ara.pred.df <- unique(data.frame('shifted_time'=ara.spline.data$shifted_time,
                                     'mean_cpm'=ara.pred, 'accession'='Col0'))
    bra.pred <- stats::predict(bra.fit)
    bra.pred.df <- unique(data.frame('shifted_time'=bra.spline.data$shifted_time,
                                     'mean_cpm'=bra.pred, 'accession'='Ro18'))

    combined.pred <- stats::predict(combined.fit)
    combined.pred.df <- unique(data.frame('shifted_time'=combined.spline.data$shifted_time,
                                          'mean_cpm'=combined.pred, 'accession'='registered'))
    spline.df <- rbind(ara.pred.df, bra.pred.df, combined.pred.df)

    p <- ggplot2::ggplot(data=combined.spline.data)+
      ggplot2::aes(x=shifted_time, y=mean_cpm, colour=accession) +
      ggplot2::geom_point()+
      ggplot2::geom_line(data=spline.df)+
      ggplot2::ggtitle(paste0(curr.sym, ' : sep AIC:combo AIC=', round(seperate.AIC), ':', round(combined.AIC),
                              ', sep BIC: combo BIC=', round(seperate.BIC), ':', round(combined.BIC)))
    # ggplot2::ggsave(paste0('./testing/fitted_splines/', curr.sym, '_', max(ara.pred.df$shifted_time), '.pdf'))
  }

  return(list(seperate.AIC, combined.AIC, seperate.BIC,combined.BIC))
}
