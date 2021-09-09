#' Calculate scores for all candidate shifts for all genes in data frame.
#'
#' `calculate_all_best_shifts` is a function to get score for all shifts value after stretching with a particular stretch factor. The extremes shift values were defined by allowing 5 overlapping points for comparison.
#'
#' @param num_shifts Number of different shifts to be considered.
#' @param mean_df Input data.
#' @param stretch_factor Current stretch factor.
#' @param do_rescale Apply "scale" to compared points for each shift if TRUE, use original mean expression data if FALSE.
#' @param shift_extreme Approximation of maximum and minimum shifts allowed.
#' @param min_num_overlapping_points Bound the extreme allowed shifts, such than at least this many timepoints are being compared for both accessions.
#' @param testing Showing a plot of the progress if TRUE, otherwise if FALSE.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_fix Accession name of data fix.
#'
#' @return
#' @export
calculate_all_best_shifts <- function(num_shifts,
                                      mean_df,
                                      stretch_factor,
                                      do_rescale,
                                      shift_extreme,
                                      min_num_overlapping_points,
                                      testing = FALSE,
                                      accession_data_to_transform = "Col0",
                                      accession_data_fix = "Ro18") {

  # Initialize vectors
  symbols <- c()
  num_points <- c()
  all_scores_list <- rep(list(0), length(unique(mean_df$locus_name)))

  # Get the extreme shifts which can be applied to the genes
  extreme_shift <- get_extreme_shifts_for_all(
    mean_df,
    stretch_factor,
    min_num_overlapping_points,
    shift_extreme,
    accession_data_to_transform,
    accession_data_fix
  )

  min_shift <- extreme_shift[[1]]
  max_shift <- extreme_shift[[2]]

  count <- 0
  for (i in 1:length(unique(mean_df$locus_name))) {
    curr_sym <- unique(mean_df$locus_name)[i]
    if (count %% 100 == 0) {
      print(paste0(count, " / ", length(unique(mean_df$locus_name))))
    }

    # Out is mean SSD between data to transform (e.g. arabidopsis), and interpolated data fix (interpolated between 2 nearest points, e.g. Brassica)
    # Get "score" for all the candidate shifts. Score is mean error / data fix expression for compared points. If timepoints don't line up, brassica value is linearly imputed
    out <- get_best_shift(num_shifts,
      curr_sym,
      data = mean_df,
      stretch_factor,
      do_rescale,
      min_shift,
      max_shift,
      testing,
      accession_data_to_transform = "Col0",
      accession_data_fix = "Ro18"
    )

    best_shift <- out$shift[out$score == min(out$score)]
    if (length(best_shift) > 1) {
      if (max(out$score) == "Inf") {
        # Can get inf score if data fix gene note expressed in the comparison
        next
      } else {
        # If ties for the best shift applied, apply the smaller absolute one
        best_shift <- best_shift[abs(best_shift) == min(abs(best_shift))]
      }
    }

    all_scores <- out
    all_scores_list[[i]] <- all_scores
    symbols <- c(symbols, curr_sym)

    count <- count + 1
  }

  all_scores_df <- do.call("rbind", all_scores_list)

  return(all_scores_df)
}



#' Calculate the score for all shifts
#'
#' `get_best_shift` is used to calculate the score for all shifts (for the current gene, and current stretch_factor), and return the scores for all as a table, and the value of the optimal shift. Shift extremes are defined s.t. at least 5 points are compared.
#'
#' @param num_shifts Number of different shifts to be considered.
#' @param curr_sym Current gene accession.
#' @param data Input data.
#' @param stretch_factor Current stretch factor.
#' @param do_rescale Apply "scale" to compared points for each shift if TRUE, use original mean expression data if FALSE.
#' @param min_shift Minimum extreme value of shift.
#' @param max_shift Maximum extreme value of shift.
#' @param testing Showing a plot of the progress if TRUE, otherwise if FALSE
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_fix Accession name of data fix.
#'
#' @export
get_best_shift <- function(num_shifts = 25,
                           curr_sym,
                           data,
                           stretch_factor,
                           do_rescale,
                           min_shift,
                           max_shift,
                           testing = FALSE,
                           accession_data_to_transform = "Col0",
                           accession_data_fix = "Ro18") {

  data <- data[data$locus_name == curr_sym, ]

  # Transform timepoint to be time from first timepoint
  data[, delta_time := timepoint - min(timepoint), by = .(accession)]

  # Apply stretch_factor to the data to transform, leave the data fix as it is
  data$delta_time[data$accession == accession_data_to_transform] <- data$delta_time[data$accession == accession_data_to_transform] * stretch_factor


  all_scores <- rep(0, num_shifts)
  all_data_transform_mean <- rep(0, num_shifts)
  all_data_fix_mean <- rep(0, num_shifts)
  all_data_transform_sd <- rep(0, num_shifts)
  all_data_fix_sd <- rep(0, num_shifts)

  all_shifts <- seq(min_shift, max_shift, length.out = num_shifts)

  if (!(0 %in% all_shifts)) {
    # Include 0 shift in candidates
    all_shifts <- c(all_shifts, 0)
  }

  # Start the iteration to calculate score for each shift for all shifts in the list
  i <- 1
  for (i in 1:length(all_shifts)) {
    curr_shift <- all_shifts[i]

    # Shift the data to transform expression timings
    data$shifted_time <- data$delta_time
    data$shifted_time[data$accession == accession_data_to_transform] <- data$delta_time[data$accession == accession_data_to_transform] + curr_shift

    # Cut down to just the data to transform and data fix timepoints which compared
    data <- get_compared_timepoints(data)
    compared <- data[data$is_compared == TRUE, ]

    # Renormalise expression using just these timepoints?
    if (do_rescale == TRUE) {
      # Record the mean and sd of the compared points, used for rescaling in "apply shift" function
      data_transform_mean <- mean(compared$mean_cpm[compared$accession == accession_data_to_transform])
      data_fix_mean <- mean(compared$mean_cpm[compared$accession == accession_data_fix])
      data_transform_sd <- stats::sd(compared$mean_cpm[compared$accession == accession_data_to_transform])
      data_fix_sd <- stats::sd(compared$mean_cpm[compared$accession == accession_data_fix])

      # Do the transformation started from here
      if ((data_transform_sd != 0 | !is.nan(data_transform_sd)) & (data_fix_sd != 0 | !is.nan(data_fix_sd))) {
        # If neither are 0, so won't be dividing by 0 (which gives NaNs)
        compared[, mean_cpm := scale(mean_cpm, scale = TRUE, center = TRUE), by = .(accession)]
      } else { # If at least one of them is all 0
        data_transform_compared <- compared[compared$accession == accession_data_to_transform, ]
        data_fix_compared <- compared[compared$accession == accession_data_fix, ]
        if ((data_transform_sd == 0) & (data_fix_sd != 0 | !is.nan(data_fix_sd))) {
          # If only data_transform_sd==0
          data_fix_compared[, mean_cpm := scale(mean_cpm, scale = TRUE, center = TRUE), by = .(accession)]
        }
        if ((data_transform_sd != 0 | !is.nan(data_transform_sd)) & (data_fix_sd == 0)) {
          # If only data_fix_sd == 0
          data_transform_compared[, mean_cpm := scale(mean_cpm, scale = TRUE, center = TRUE), by = .(accession)]
        }
        # If both are all 0, then do nothing.
        compared <- rbind(data_transform_compared, data_fix_compared)
      }
    } else {
      # If didn't rescale expression for comparison, record values s.t. (x - xmean) / x_sd = x
      data_transform_mean <- 0
      data_fix_mean <- 0
      data_transform_sd <- 1
      data_fix_sd <- 1
    }

    # Data plot of shifted, and normalised gene expression
    # if (testing == TRUE) {
    #   p <- ggplot2::ggplot(compared) +
    #     ggplot2::aes(x = shifted_time, y = mean_cpm, color = accession) +
    #     ggplot2::geom_point() +
    #     ggplot2::geom_line() +
    #     ggplot2::ggtitle(paste0("shift : ", curr_shift))
    #
    #   ggplot2::ggsave(
    #     plot = p,
    #     filename = paste0(curr_sym, stretch_factor, "-", curr_shift, ".pdf")
    #   )
    # }

    # For each data to transform timepoint, linear interpolate between the two nearest data fix timepoints
    data_transform_compared <- compared[compared$accession == accession_data_to_transform]
    data_fix_compared <- compared[compared$accession == accession_data_fix]

    data_transform_compared$pred.bra.expression <- sapply(data_transform_compared$shifted_time, interpolate_data_fix_comparison_expression, data_fix_dt = data_fix_compared)

    # if (testing == TRUE) {
    #   interpolate_res <- ggplot2::ggplot(compared) +
    #     ggplot2::aes(x = shifted_time, y = mean_cpm, color = accession) +
    #     ggplot2::geom_point() +
    #     ggplot2::geom_line() +
    #     ggplot2::geom_point(data = data_transform_compared, aes(x = shifted_time, y = pred.bra.expression), color = "purple") +
    #     ggplot2::geom_line(data = data_transform_compared, aes(x = shifted_time, y = pred.bra.expression), color = "purple")
    #
    #
    #     ggplot2::ggsave(
    #       plot = interpolate_res,
    #       filename = paste0(curr_sym, stretch_factor, "-", curr_shift, "with_interpolation.pdf")
    #     )
    # }

    # Calculate the score, using the (interpolated) predicted.bra.expression, and the observed arabidopsis expression
    score <- calc_score(data_transform_compared$mean_cpm, data_transform_compared$pred.bra.expression)

    if (is.na(score)) {
      print("error in get_best_shift(): got a score of NA for gene:")
      print(paste("data_transform_compared$mean_cpm:", data_transform_compared$mean_cpm))
      print(data_transform_compared$pred.bra.expression)
      print(curr_sym)
      print(paste0("with curr_shift=", curr_shift))
      stop()
    }

    all_scores[i] <- score
    all_data_transform_mean[i] <- data_transform_mean
    all_data_fix_mean[i] <- data_fix_mean
    all_data_transform_sd[i] <- data_transform_sd
    all_data_fix_sd[i] <- data_fix_sd
  }

  out <- data.table::data.table(data.frame(
    "gene" = curr_sym,
    "stretch" = stretch_factor,
    "shift" = all_shifts,
    "score" = all_scores,
    "data_transform_compared_mean" = all_data_transform_mean,
    "data_fix_compared_mean" = all_data_fix_mean,
    "data_transform_compared_sd" = all_data_transform_sd,
    "data_fix_compared_sd" = all_data_fix_sd
  ))

  return(out)

}
