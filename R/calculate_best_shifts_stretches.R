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
#' @param accession_data_to_align Accession name of data which will be aligned.
#' @param accession_data_target Accession name of data target.
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
                                      accession_data_to_align = "Col0",
                                      accession_data_target = "Ro18") {

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
    accession_data_to_align,
    accession_data_target
  )

  min_shift <- extreme_shift[[1]]
  max_shift <- extreme_shift[[2]]

  count <- 0
  for (i in 1:length(unique(mean_df$locus_name))) {
    curr_sym <- unique(mean_df$locus_name)[i]
    if (count %% 100 == 0) {
      print(paste0(count, " / ", length(unique(mean_df$locus_name))))
    }

    # Out is mean SSD between data to align (e.g. arabidopsis), and interpolated data target (interpolated between 2 nearest points, e.g. Brassica)
    # Get "score" for all the candidate shifts. Score is mean error / data target expression for compared points. If timepoints don't line up, brassica value is linearly imputed
    out <- get_best_shift(num_shifts,
      curr_sym,
      data = mean_df,
      stretch_factor,
      do_rescale,
      min_shift,
      max_shift,
      testing = FALSE,
      accession_data_to_align = "Col0",
      accession_data_target = "Ro18"
    )

    best_shift <- out$shift[out$score == min(out$score)]
    if (length(best_shift) > 1) {
      if (max(out$score) == "Inf") {
        # Can get inf score if data target gene note expressed in the comparison
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
#' @param accession_data_to_align Accession name of data which will be aligned.
#' @param accession_data_target Accession name of data target.
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
                           accession_data_to_align = "Col0",
                           accession_data_target = "Ro18") {

  data <- data[data$locus_name == curr_sym, ]

  # Transform timepoint to be time from first timepoint
  data[, delta_time := timepoint - min(timepoint), by = .(accession)]

  # Apply stretch_factor to the data to align, leave the data target as it is
  data$delta_time[data$accession == accession_data_to_align] <- data$delta_time[data$accession == accession_data_to_align] * stretch_factor


  all_scores <- rep(0, num_shifts)
  all_data_align_mean <- rep(0, num_shifts)
  all_data_target_mean <- rep(0, num_shifts)
  all_data_align_sd <- rep(0, num_shifts)
  all_data_target_sd <- rep(0, num_shifts)

  all_shifts <- seq(min_shift, max_shift, length.out = num_shifts)

  if (!(0 %in% all_shifts)) {
    # Include 0 shift in candidates
    all_shifts <- c(all_shifts, 0)
  }

  # Start the iteration to calculate score for each shift for all shifts in the list
  i <- 1
  for (i in 1:length(all_shifts)) {
    curr_shift <- all_shifts[i]

    # Shift the data to align expression timings
    data$shifted_time <- data$delta_time
    data$shifted_time[data$accession == accession_data_to_align] <- data$delta_time[data$accession == accession_data_to_align] + curr_shift

    # Cut down to just the data to align and data target timepoints which compared
    data <- get_compared_timepoints(data)
    compared <- data[data$is_compared == TRUE, ]

    # Renormalise expression using just these timepoints?
    if (do_rescale == TRUE) {
      # Record the mean and sd of the compared points, used for rescaling in "apply shift" function
      data_align_mean <- mean(compared$mean_cpm[compared$accession == accession_data_to_align])
      data_target_mean <- mean(compared$mean_cpm[compared$accession == accession_data_target])
      data_align_sd <- stats::sd(compared$mean_cpm[compared$accession == accession_data_to_align])
      data_target_sd <- stats::sd(compared$mean_cpm[compared$accession == accession_data_target])

      # Do the transformation started from here
      if ((data_align_sd != 0 | !is.nan(data_align_sd)) & (data_target_sd != 0 | !is.nan(data_target_sd))) {
        # If neither are 0, so won't be dividing by 0 (which gives NaNs)
        compared[, mean_cpm := scale(mean_cpm, scale = TRUE, center = TRUE), by = .(accession)]
      } else { # If at least one of them is all 0
        data_align_compared <- compared[compared$accession == accession_data_to_align, ]
        data_target_compared <- compared[compared$accession == accession_data_target, ]
        if ((data_align_sd == 0) & (data_target_sd != 0 | !is.nan(data_target_sd))) {
          # If only data_align_sd==0
          data_target_compared[, mean_cpm := scale(mean_cpm, scale = TRUE, center = TRUE), by = .(accession)]
        }
        if ((data_align_sd != 0 | !is.nan(data_align_sd)) & (data_target_sd == 0)) {
          # If only data_target_sd == 0
          data_align_compared[, mean_cpm := scale(mean_cpm, scale = TRUE, center = TRUE), by = .(accession)]
        }
        # If both are all 0, then do nothing.
        compared <- rbind(data_align_compared, data_target_compared)
      }
    } else {
      # If didn't rescale expression for comparison, record values s.t. (x - xmean) / x_sd = x
      data_align_mean <- 0
      data_target_mean <- 0
      data_align_sd <- 1
      data_target_sd <- 1
    }

    # Data plot of shifted, and normalised gene expression
    if (testing == TRUE) {
      p <- ggplot2::ggplot(compared) +
        ggplot2::aes(x = shifted_time, y = mean_cpm, color = accession) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::ggtitle(paste0("shift : ", curr_shift))
      ggplot2::ggsave(paste0(curr_sym, stretch_factor, "-", curr_shift, ".pdf"))
    }

    # For each data to align timepoint, linear interpolate between the two nearest data target timepoints
    data_align_compared <- compared[compared$accession == accession_data_to_align]
    data_target_compared <- compared[compared$accession == accession_data_target]

    data_align_compared$pred.bra.expression <- sapply(data_align_compared$shifted_time, interpolate_data_target_comparison_expression, data_target_dt = data_target_compared)

    if (testing == TRUE) {
      interpolate_res <- ggplot2::ggplot(compared) +
        ggplot2::aes(x = shifted_time, y = mean_cpm, color = accession) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::geom_point(data = data_align_compared, aes(x = shifted_time, y = pred.bra.expression), color = "purple") +
        ggplot2::geom_line(data = data_align_compared, aes(x = shifted_time, y = pred.bra.expression), color = "purple") +
        ggplot2::ggsave(paste0(curr_sym, stretch_factor, "-", curr_shift, "with_interpolation.pdf"))
    }

    # Calculate the score, using the (interpolated) predicted.bra.expression, and the observed arabidopsis expression
    score <- calc_score(data_align_compared$mean_cpm, data_align_compared$pred.bra.expression)

    if (is.na(score)) {
      print("error in get_best_shift(): got a score of NA for gene:")
      print(paste("data_align_compared$mean_cpm:", data_align_compared$mean_cpm))
      print(data_align_compared$pred.bra.expression)
      print(curr_sym)
      print(paste0("with curr_shift=", curr_shift))
      stop()
    }

    all_scores[i] <- score
    all_data_align_mean[i] <- data_align_mean
    all_data_target_mean[i] <- data_target_mean
    all_data_align_sd[i] <- data_align_sd
    all_data_target_sd[i] <- data_target_sd
  }

  out <- data.table::data.table(data.frame(
    "gene" = curr_sym,
    "stretch" = stretch_factor,
    "shift" = all_shifts,
    "score" = all_scores,
    "data_align_compared_mean" = all_data_align_mean,
    "data_target_compared_mean" = all_data_target_mean,
    "data_align_compared_sd" = all_data_align_sd,
    "data_target_compared_sd" = all_data_target_sd
  ))

  return(out)

}
