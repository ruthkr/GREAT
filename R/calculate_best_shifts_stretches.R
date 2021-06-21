# for all genes, in mean_df, get the scores, and the shifts which result in them after
# stretching arabidopsis by "stretch_factor, and applying shifts forward and
# backward, with the extremes defined by allowing 5 overlapping points for comparison.
# mean_df <- to.shift.df
# stretch_factor <- 2
# do_rescale=F
#' @export
calculate_all_best_shifts <- function(mean_df,
                                      stretch_factor,
                                      min_num_overlapping_points,
                                      shift_extreme,
                                      accession_data_to_align,
                                      accession_data_target,
                                      do_rescale) {

  # Initialize vectors
  symbols <- c()
  num_points <- c()
  all_scores_list <- rep(list(0), length(unique(mean_df$locus_name)))

  # Get the extreme shifts which can be applied to the genes
  extreme_shift <- get_extreme_shifts_for_all(mean_df,
    stretch_factor,
    min_num_overlapping_points,
    shift_extreme,
    accession_data_to_align,
    accession_data_target)

  min_shift <- extreme_shift[[1]]
  max_shift <- extreme_shift[[2]]

  print(paste0("min shift:", min_shift, "max shift:", max_shift))

  count <- 0
  for (i in 1:length(unique(mean_df$locus_name))) {
    curr_sym <- unique(mean_df$locus_name)[i]
    if (count %% 100 == 0) {
      print(paste0(count, ' / ', length(unique(mean_df$locus_name))))
    }

    # out is mean SSD between arabidopsis, and interpolated brassica (interpolated between 2 nearest points)
    # ggplot2::ggplot(mean_df[mean_df$locus_name==curr_sym,])+
    #   ggplot2::aes(x=timepoint, y=mean_cpm, color=accession) +
    #   ggplot2::geom_point()

    ### get "score" for all the candidate shifts - score is mean error / brassica expression for compared points.
    ### if timepoints don't line up, brassica value is linearly imputed
    out <- get_best_shift(curr_sym, mean_df, stretch_factor, do_rescale, min_shift, max_shift, testing=FALSE)

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



# wroking on new implementation of the funciton
# curr_sym <- 'BRAA02G015410.3C'
#test <- mean_df
# stretch_factor <- 1
# do_rescale <- TRUE
#testing=T

#' @export
get_best_shift <- function(curr_sym,
                           test,
                           stretch_factor,
                           do_rescale,
                           min_shift,
                           max_shift,
                           testing = FALSE,
                           accession_data_to_align = "Col0",
                           accession_data_target = "Ro18") {

  # for the current gene, and current stretch_factor, calculate the score for all
  # shifts, and return the scores for all as a table, and the value of the optimal shift.

  # Shift extremes are defined s.t. at least 5 points are compared.

  # do_rescale == TRUE, means apply "scale" to compared points for each shift. ==FALSE, means use original mean expression data

  num_shifts <- 25 # the number of different shifts to be considered.


  test <- test[test$locus_name==curr_sym, ]

  # transform timepoint to be time from first timepoint
  test[, delta_time:=timepoint - min(timepoint), by=.(accession)]
  # apply stretch_factor to the arabidopsis, leave the rapa as is
  test$delta_time[test$accession == 'Col0'] <- test$delta_time[test$accession=='Col0'] * stretch_factor


  # # Try to make it consistent as in apply stretch
  # test$delta_time[test$accession=='Col0'] <- test$delta_time[test$accession=='Col0'] + 7
  # test$delta_time[test$accession=='Ro18'] <- test$delta_time[test$accession=='Ro18'] + 14

  all_scores <- rep(0, num_shifts)
  all.ara.mean <- rep(0, num_shifts)
  all.bra.mean <- rep(0, num_shifts)
  all.ara.sd <- rep(0, num_shifts)
  all.bra.sd <- rep(0, num_shifts)

  all_shifts <- seq(min_shift, max_shift, length_out=num_shifts)
  if (!(0 %in% all_shifts)) {
    all_shifts <- c(all_shifts, 0) # include 0 shift in candidates.
  }
  i=1
  for (i in 1:length(all_shifts)) {

    curr_shift <- all_shifts[i]

    # shift the arabidopsis expression timings
    test$shifted_time <- test$delta_time
    test$shifted_time[test$accession == 'Col0'] <- test$delta_time[test$accession == 'Col0'] + curr_shift

    #### test plot - of shifted, UNNORMALISED gene expression
    # if (testing==TRUE) {
    #   p <- ggplot2::ggplot(test)+
    #     ggplot2::aes(x=shifted_time, y=mean_cpm, color=accession)+
    #     ggplot2::geom_point()+
    #     ggplot2::ggtitle(paste0('shift : ', curr_shift))
    #   p
    #   ggplot2::ggsave(paste0('./testing/', curr_shift, '.pdf'))
    # }


    # Cut down to just the arabidopsis and brassica timepoints which compared
    test <- get_compared_timepoints(test)
    compared <- test[test$is.compared==TRUE, ]

    # Renormalise expression using just these timepoints?
    if (do_rescale == TRUE) {
      # record the mean and sd of the compared points, used for rescaling
      # in "apply shift" function
      ara.mean <- mean(compared$mean_cpm[compared$accession=='Col0'])
      bra.mean <- mean(compared$mean_cpm[compared$accession=='Ro18'])
      ara.sd <- stats::sd(compared$mean_cpm[compared$accession=='Col0'])
      bra.sd<- stats::sd(compared$mean_cpm[compared$accession=='Ro18'])

      # do the transformation for here
      if ((ara.sd != 0 | !is.nan(ara.sd)) & (bra.sd != 0 | !is.nan(bra.sd))) { # if neither are 0, so won't be dividing by 0 (which gives NaNs)
        compared[, mean_cpm:=scale(mean_cpm, scale=TRUE, center=TRUE), by=.(accession)]
      } else { # if at least one of them is all 0
        ara.compared <- compared[compared$accession=='Col0',]
        bra.compared <- compared[compared$accession=='Ro18',]
        if((ara.sd == 0) & (bra.sd != 0 | !is.nan(bra.sd))) { # if only ara.sd==0
          bra.compared[, mean_cpm:=scale(mean_cpm, scale=TRUE, center=TRUE), by=.(accession)]
        }
        if ((ara.sd != 0 | !is.nan(ara.sd))  & (bra.sd == 0)) { # if only bra.sd == 0
          ara.compared[, mean_cpm:=scale(mean_cpm, scale=TRUE, center=TRUE), by=.(accession)]
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
        ggplot2::aes(x=shifted_time, y=mean_cpm, color=accession)+
        ggplot2::geom_point()+
        ggplot2::ggtitle(paste0('shift : ', curr_shift))
      ggplot2::ggsave(paste0('./testing/',stretch_factor, '-', curr_shift, '.pdf'))
    }

    # for each arabidopsis timepoint, linear interpolate between the two nearest brassica timepoints
    ara.compared <- compared[compared$accession=='Col0']
    bra.compared <- compared[compared$accession=='Ro18']

    ara.compared$pred.bra.expression <- sapply(ara.compared$shifted_time, interpolate_brassica_comparison_expression, bra.dt=bra.compared)

    # calculate the score, using the (interpolated) predicted.bra.expression, and the observed arabidopsis expression
    # score = mean ((observed - expected)**2 )
    score <- calc_score(ara.compared$mean_cpm, ara.compared$pred.bra.expression)

    if (is.na(score)) {
      print('error in get_best_shift(): got a score of NA for gene:')
      print(paste("ara.compared$mean_cpm:", ara.compared$mean_cpm))
      print(ara.compared$pred.bra.expression)
      print(curr_sym)
      print(paste0('with curr_shift=', curr_shift))
      stop()
    }

    all_scores[i] <- score
    all.ara.mean[i] <- ara.mean
    all.bra.mean[i] <- bra.mean
    all.ara.sd[i] <- ara.sd
    all.bra.sd[i] <- bra.sd
  }

  out <- data.table::data.table(data.frame('gene'=curr_sym, 'stretch'=stretch_factor, 'shift'=all_shifts, 'score'=all_scores,
                                           'ara.compared.mean'=all.ara.mean, 'bra.compared.mean'=all.bra.mean,
                                           'ara.compared.sd'=all.ara.sd, 'bra.compared.sd'=all.bra.sd))
  return(out)
}
