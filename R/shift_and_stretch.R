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


#' Wrapper of applying best shifts and compare the registered and unregistered models
#'
#' `calculate_all_model_comparison_stats` is a wrapper function to apply best shifts and compare the registered and unregistered models using compare_registered_to_unregistered_model.
#'
#' @param all_data_df Input all data (without taking mean).
#' @param best_shifts Input data frame containing information of best shifts.
#' @param accession_data_to_align Accession name of data which will be aligned.
#' @param accession_data_target Accession name of data target.
#' @param data_to_align_time_added Time points to be added in data to align.
#' @param data_target_time_added Time points to be added in data target.
#'
#' @return AIC and BIC score for registered and unregistered models.
#' @export
calculate_all_model_comparison_stats <- function(all_data_df,
                                                 best_shifts,
                                                 accession_data_to_align,
                                                 accession_data_target,
                                                 data_to_align_time_added,
                                                 data_target_time_added) {

  if (!(accession_data_to_align %in% unique(all_data_df$accession) &
        accession_data_target %in% unique(all_data_df$accession))) {
    stop("error in calculate_all_model_comparison_stats() :
         all_data_df doesn't have the correct accession info - should have been
         converted to Ro18 & Col0")
  }


  # Apply the registration to the all rep data, so can use for model comparison
  shifted_all_data_df <- apply_best_shift(data = all_data_df,
                                          best_shifts,
                                          accession_data_to_align,
                                          accession_data_target,
                                          data_to_align_time_added,
                                          data_target_time_added)

  message('Calculating registration vs different expression comparison AIC & BIC...')

  genes <- unique(shifted_all_data_df$locus_name)

  out.sepAIC <- rep(0, length(genes))
  out.combAIC <- rep(0, length(genes))
  out.sepBIC <- rep(0, length(genes))
  out.combBIC <- rep(0, length(genes))

  for(i in 1:length(genes)) {

    if( i %% 100 == 0) {
      message(i, ' / ', length(genes))
    }

    curr_sym <- genes[i]

    L <- compare_registered_to_unregistered_model(curr_sym,
                                                  shifted_all_data_df,
                                                  is_testing = FALSE,
                                                  accession_data_to_align,
                                                  accession_data_target)

    out.sepAIC[i] <- L[["seperate.AIC"]]
    out.combAIC[i] <- L[["combined.AIC"]]
    out.sepBIC[i] <- L[["seperate.BIC"]]
    out.combBIC[i] <- L[["combined.BIC"]]

  }

  out <- data.table::data.table('gene' = genes,
                                           'seperate.AIC' = out.sepAIC,
                                           'registered.AIC' = out.combAIC,
                                           'seperate.BIC' = out.sepBIC,
                                          'registered.BIC' = out.combBIC)

  return(out)

}


#' Register all expression over time using optimal shift found
#'
#' `apply_best_shift` is a function to register all unregistered expression overtime using the best/optimal shifts found.
#'
#' @param data Input data (all data).
#' @param best_shifts Input data frame containing information of best shifts.
#' @param accession_data_to_align Accession name of data which will be aligned.
#' @param accession_data_target Accession name of data target.
#' @param data_to_align_time_added Time points to be added in data to align.
#' @param data_target_time_added Time points to be added in data target.
#'
#' @return The registered expression over time for each gene.
#' @export
apply_best_shift <- function(data,
                             best_shifts,
                             accession_data_to_align,
                             accession_data_target,
                             data_to_align_time_added,
                             data_target_time_added) {

  processed_data <- data.table::copy(data)

  processed_data <- apply_stretch(
    data = processed_data,
    best_shifts,
    accession_data_to_align,
    accession_data_target,
    data_to_align_time_added,
    data_target_time_added
  )

  # Normalise the expression data (If was normalised when calculating the expression data, is recorder in the _compared_mean, and _compared_sd columns. If no normalisation was carried out, then these should have values of 0 and 1. This was done using get_best_shift()).


  if (!(all(unique(best_shifts$data_align_compared_mean) == 0)) | !(all(unique(best_shifts$data_target_compared_mean) == 0))) {

    message("Normalising expression by mean and sd of compared values...")

    processed_data <- apply_best_normalisation(
      data = processed_data,
      best_shifts,
      accession_data_to_align,
      accession_data_target
    )

    message("Done!")

  } else {

    # If no scaling carried out DURING the registration step
    message("No normalisation was carried out DURING registration (though may have been, prior to Col-Ro18 comparison)")

    processed_data <- processed_data

  }

  message("Applying best shift...")

  # For each gene, shift the data to align expression by the optimal shift found previously

  for (curr_gene in unique(processed_data$locus_name)) {

    curr_best_shift <- best_shifts$shift[best_shifts$gene == curr_gene]
    processed_data$shifted_time[processed_data$accession == accession_data_to_align & processed_data$locus_name == curr_gene] <- processed_data$shifted_time[processed_data$accession == accession_data_to_align & processed_data$locus_name == curr_gene] + curr_best_shift

  }

  message("Done!")

  return(processed_data)

}


#' Apply stretch factor
#'
#' @param data Input data.
#' @param best_shifts Input data frame containing information of best shifts.
#' @param accession_data_to_align Accession name of data which will be aligned.
#' @param accession_data_target Accession name of data target.
#' @param data_to_align_time_added Time points to be added in data to align.
#' @param data_target_time_added Time points to be added in data target.
#'
#' @return
#' @export
apply_stretch <- function(data,
                          best_shifts,
                          accession_data_to_align = "Col0",
                          accession_data_target = "Ro18",
                          data_to_align_time_added = 11,
                          data_target_time_added = 11) {

  data <- data.table::copy(data)

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



#' Apply normalisation (after applying stretch)
#'
#' `apply_best_normalisation` is a function to normalise by the mean and standard deviation of the compared points (after applying stretching) for each gene, in each accesion (data target and data to align). If the gene wasn't compared, set the expression value to NA.
#'
#' @param data Input data (after applying stretching).
#' @param best_shifts Input dataframe containing information of best shifts.
#' @param accession_data_to_align Accession name of data which will be aligned.
#' @param accession_data_target Accession name of data target.
#'
#' @return Normalised data.
#' @export
apply_best_normalisation <- function(data,
                                     best_shifts,
                                     accession_data_to_align = "Col0",
                                     accession_data_target = "Ro18") {

  count <- 0
  for (curr_gene in unique(data$locus_name)) {

    if (count %% 100 == 0) {
      message(count, " / ", length(unique(data$locus_name)))
    }

    data_align_mean <- best_shifts$data_align_compared_mean[best_shifts$gene == curr_gene]
    data_target_mean <- best_shifts$data_target_compared_mean[best_shifts$gene == curr_gene]
    data_align_sd <- best_shifts$data_align_compared_sd[best_shifts$gene == curr_gene]
    data_target_sd <- best_shifts$data_target_compared_sd[best_shifts$gene == curr_gene]

    # If was compared
    if (length(data_align_mean) != 0) {

      # Make sure that sd is not 0, since we do not want to divide by 0 ---------------
      if (data_align_sd != 0) {
        data$mean_cpm[data$locus_name == curr_gene & data$accession == accession_data_to_align] <- (data$mean_cpm[data$locus_name == curr_gene & data$accession == accession_data_to_align] - data_align_mean) / data_align_sd
      } else {
        data$mean_cpm[data$locus_name == curr_gene & data$accession == accession_data_to_align] <- (data$mean_cpm[data$locus_name == curr_gene & data$accession == accession_data_to_align] - data_align_mean)
      }

      # Make sure that sd is not 0, since we do not want to divide by 0 ---------------
      if (data_target_sd != 0) {
        data$mean_cpm[data$locus_name == curr_gene & data$accession == accession_data_target] <- (data$mean_cpm[data$locus_name == curr_gene & data$accession == accession_data_target] - data_target_mean) / data_target_sd
      } else {
        data$mean_cpm[data$locus_name == curr_gene & data$accession == accession_data_target] <- (data$mean_cpm[data$locus_name == curr_gene & data$accession == accession_data_target] - data_target_mean)
      }

      if (any(is.na(data$mean_cpm))) {
        message("Have NAs in mean_cpm after rescaling in apply best_normalisation() for gene :")
        message(unique(data$locus_name))
        stop()
      }

    } else {

      data$mean_cpm[data$locus_name == curr_gene & data$accession == accession_data_to_align] <- NA
      data$mean_cpm[data$locus_name == curr_gene & data$accession == accession_data_target] <- NA

    }

    count <- count + 1

  }

  return(data)

}



#' Comparing registered to unregistered model
#'
#' `compare_registered_to_unregistered_model` is a function to compare the overlapping timepoints in data target and data to align after best registration and without registration. Same timepoints and stretched data were used for both models.
#'
#' @param curr_sym A gene accession.
#' @param all_data_df Input data.
#' @param is_testing Showing a plot of the progress if TRUE, otherwise if FALSE.
#' @param accession_data_to_align Accession name of data which will be aligned.
#' @param accession_data_target Accession name of data target.
#'
#' @return Score of AIC and BIC for both registered and unregistered models.
#' @export
compare_registered_to_unregistered_model <- function(curr_sym,
                                                     all_data_df,
                                                     is_testing = FALSE,
                                                     accession_data_to_align = "Col0",
                                                     accession_data_target = "Ro18") {


  curr_data_df <- all_data_df[all_data_df$locus_name == curr_sym]


  # Flag the timepoints to be used in the modelling, only the ones which overlap!
  curr_data_df <- get_compared_timepoints(curr_data_df,
                                          accession_data_to_align,
                                          accession_data_target)


#   ggplot2::ggplot(curr_data_df)+
#     ggplot2::aes(x=shifted_time, y=mean_cpm, shape=is_compared, color=accession)+
#     ggplot2::geom_point()

  # Cut down to the data for each model
  data_to_align_spline <- curr_data_df[curr_data_df$is_compared == TRUE &
                                    curr_data_df$accession == accession_data_to_align, ]
  data_target_spline <- curr_data_df[curr_data_df$is_compared == TRUE &
                                    curr_data_df$accession == accession_data_target, ]
  combined_spline_data <- curr_data_df[curr_data_df$is_compared == TRUE, ]

  # Fit the models - fit regression splines.
  # http://www.utstat.utoronto.ca/reid/sta450/feb23.pdf
  # for cubic spline, K+3 params where K=num.knots
  # as can omit constant term
  num.spline.params <- 6 # number of parameters for each spline fitting (degree and this used to calculate num knots).
  num.registration.params <- 2 # stretch, shift
  num.obs <- nrow(combined_spline_data)


  data_to_align_fit <- stats::lm(mean_cpm ~ splines::bs(shifted_time, df = num.spline.params, degree = 3), data = data_to_align_spline)
  data_target_fit <- stats::lm(mean_cpm ~ splines::bs(shifted_time, df = num.spline.params, degree = 3), data = data_target_spline)
  combined_fit <- stats::lm(mean_cpm ~ splines::bs(shifted_time, df = num.spline.params, degree = 3), data = combined_spline_data)

  # Calculate the log likelihoods
  data_to_align_logLik <- stats::logLik(data_to_align_fit)
  data_target_logLik <- stats::logLik(data_target_fit)
  seperate_logLik <- data_to_align_logLik + data_target_logLik # logLikelihoods, so sum
  combined_logLik <- stats::logLik(combined_fit)

  # Calculate the comparison.stats - - AIC, BIC, smaller is better!
  # 2*num.spline.params as fitting seperate models for Ara * Col
  seperate.AIC <- calc_AIC(seperate_logLik, 2*num.spline.params)
  combined.AIC <- calc_AIC(combined_logLik, num.spline.params+num.registration.params)

  seperate.BIC <- calc_BIC(seperate_logLik, 2*num.spline.params, num.obs)
  combined.BIC <- calc_BIC(combined_logLik, num.spline.params+num.registration.params, num.obs)


  if (is_testing == TRUE) {
    ara.pred <- stats::predict(data_to_align_fit)
    ara.pred.df <- unique(data.frame('shifted_time'=data_to_align_spline$shifted_time,
                                     'mean_cpm'=ara.pred, 'accession'='Col0'))
    bra.pred <- stats::predict(data_target_fit)
    bra.pred.df <- unique(data.frame('shifted_time'=data_target_spline$shifted_time,
                                     'mean_cpm'=bra.pred, 'accession'='Ro18'))

    combined.pred <- stats::predict(combined_fit)
    combined.pred.df <- unique(data.frame('shifted_time'=combined_spline_data$shifted_time,
                                          'mean_cpm'=combined.pred, 'accession'='registered'))
    spline.df <- rbind(ara.pred.df, bra.pred.df, combined.pred.df)

    p <- ggplot2::ggplot(data=combined_spline_data)+
      ggplot2::aes(x=shifted_time, y=mean_cpm, colour=accession) +
      ggplot2::geom_point()+
      ggplot2::geom_line(data=spline.df)+
      ggplot2::ggtitle(paste0(curr_sym, ' : sep AIC:combo AIC=', round(seperate.AIC), ':', round(combined.AIC),
                              ', sep BIC: combo BIC=', round(seperate.BIC), ':', round(combined.BIC)))
    ggplot2::ggsave(paste0(curr_sym, '_', max(ara.pred.df$shifted_time), '.pdf'))

  }

  return(list(seperate.AIC = seperate.AIC,
              combined.AIC = combined.AIC,
              seperate.BIC = seperate.BIC,
              combined.BIC = combined.BIC))

}
