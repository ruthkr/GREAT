#' Calculate scores for all candidate shifts for all genes in data frame.
#'
#' @noRd
calculate_all_best_shifts <- function(num_shifts,
                                      mean_df,
                                      stretch_factor,
                                      do_rescale,
                                      shift_extreme,
                                      min_num_overlapping_points,
                                      accession_data_to_transform,
                                      accession_data_ref) {
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
    accession_data_ref
  )

  min_shift <- extreme_shift[[1]]
  max_shift <- extreme_shift[[2]]

  i <- 0
  cli::cli_progress_step("Calculating score for all shifts ({i}/{length(unique(mean_df$locus_name))})", spinner = TRUE)
  for (i in seq_along(unique(mean_df$locus_name))) {
    curr_sym <- unique(mean_df$locus_name)[i]

    # Out is mean SSD between data to transform (e.g. arabidopsis), and interpolated reference data (interpolated between 2 nearest points, e.g. Brassica)
    # Get "score" for all the candidate shifts. Score is mean error / reference data expression for compared points. If timepoints don't line up, brassica value is linearly imputed
    out <- get_best_shift(
      num_shifts,
      curr_sym,
      data = mean_df,
      stretch_factor,
      do_rescale,
      min_shift,
      max_shift,
      accession_data_to_transform,
      accession_data_ref
    )

    best_shift <- out$shift[out$score == min(out$score)]

    if (length(best_shift) > 1) {
      if (is.infinite(max(out$score))) {
        # Can get inf score if reference data gene not expressed in the comparison
        next
      } else {
        # If ties for the best shift applied, apply the smaller absolute one
        best_shift <- best_shift[abs(best_shift) == min(abs(best_shift))]
      }
    }

    all_scores <- out
    all_scores_list[[i]] <- all_scores
    symbols <- c(symbols, curr_sym)

    cli::cli_progress_update()
  }

  # Bind all scores
  all_scores_df <- do.call("rbind", all_scores_list)

  return(all_scores_df)
}

#' Calculate the score for all shifts
#'
#' @noRd
get_best_shift <- function(num_shifts = 25,
                           curr_sym,
                           data,
                           stretch_factor,
                           do_rescale,
                           min_shift,
                           max_shift,
                           accession_data_to_transform,
                           accession_data_ref) {
  # Suppress "no visible binding for global variable" note
  delta_time <- NULL
  timepoint <- NULL
  accession <- NULL
  expression_value <- NULL

  # Filter locus_name
  data <- data[data$locus_name == curr_sym, ]

  # Transform timepoint to be time from first timepoint
  data[, delta_time := timepoint - min(timepoint), by = .(accession)]

  # Apply stretch_factor to the data to transform, leave the reference data as it is
  data$delta_time[data$accession == accession_data_to_transform] <- data$delta_time[data$accession == accession_data_to_transform] * stretch_factor

  all_scores <- numeric(length = num_shifts)
  all_data_transform_mean <- numeric(length = num_shifts)
  all_data_ref_mean <- numeric(length = num_shifts)
  all_data_transform_sd <- numeric(length = num_shifts)
  all_data_ref_sd <- numeric(length = num_shifts)

  all_shifts <- seq(min_shift, max_shift, length.out = num_shifts)

  if (!(0 %in% all_shifts)) {
    # Include 0 shift in candidates
    all_shifts <- c(all_shifts, 0)
  }

  # Start the iteration to calculate score for each shift for all shifts in the list
  i <- 1
  for (i in seq_along(all_shifts)) {
    curr_shift <- all_shifts[i]

    # Shift the data to transform expression timings
    data$shifted_time <- data$delta_time
    data$shifted_time[data$accession == accession_data_to_transform] <- data$delta_time[data$accession == accession_data_to_transform] + curr_shift

    # Cut down to just the data to transform and reference data timepoints which compared
    data <- get_compared_timepoints(
      data,
      accession_data_to_transform,
      accession_data_ref
    )
    compared <- data[data$is_compared == TRUE, ]

    # Renormalise expression using just these timepoints?
    if (do_rescale) {
      # Record the mean and sd of the compared points, used for rescaling in "apply shift" function
      expression_value_transform <- compared$expression_value[compared$accession == accession_data_to_transform]
      data_transform_mean <- mean(expression_value_transform)
      data_transform_sd <- stats::sd(expression_value_transform)
      expression_value_ref <- compared$expression_value[compared$accession == accession_data_ref]
      data_ref_mean <- mean(expression_value_ref)
      data_ref_sd <- stats::sd(expression_value_ref)

      # Do the transformation started from here
      if ((data_transform_sd != 0 & !is.nan(data_transform_sd)) & (data_ref_sd != 0 & !is.nan(data_ref_sd))) {
        # If neither are 0, so won't be dividing by 0 (which gives NaNs)
        compared[, expression_value := scale(expression_value, scale = TRUE, center = TRUE), by = .(accession)]
      } else { # If at least one of them is all 0
        data_transform_compared <- compared[compared$accession == accession_data_to_transform, ]
        data_ref_compared <- compared[compared$accession == accession_data_ref, ]
        if ((data_transform_sd == 0) & (data_ref_sd != 0 & !is.nan(data_ref_sd))) {
          # If only data_transform_sd==0
          data_ref_compared[, expression_value := scale(expression_value, scale = TRUE, center = TRUE), by = .(accession)]
        }
        if ((data_transform_sd != 0 & !is.nan(data_transform_sd)) & (data_ref_sd == 0)) {
          # If only data_ref_sd == 0
          data_transform_compared[, expression_value := scale(expression_value, scale = TRUE, center = TRUE), by = .(accession)]
        }
        # If both are all 0, then do nothing.
        compared <- rbind(data_transform_compared, data_ref_compared)
      }
    } else {
      cli::cli_alert_info("Not rescaling in get_best_shift()!")
      # If didn't rescale expression for comparison, record values s.t. (x - xmean) / x_sd = x
      data_transform_mean <- 0
      data_ref_mean <- 0
      data_transform_sd <- 1
      data_ref_sd <- 1
    }

    # For each data to transform timepoint, linear interpolate between the two nearest reference data timepoints
    data_transform_compared <- compared[compared$accession == accession_data_to_transform]
    data_ref_compared <- compared[compared$accession == accession_data_ref]

    data_transform_compared$pred_data_ref_expression <- sapply(data_transform_compared$shifted_time, interpolate_data_ref_comparison_expression, data_ref_dt = data_ref_compared)

    # Calculate the score, using the (interpolated) predicted.bra.expression, and the observed arabidopsis expression
    score <- calc_score(
      data_to_transform_expression = data_transform_compared$expression_value,
      data_ref_expression = data_transform_compared$pred_data_ref_expression
    )

    if (is.na(score)) {
      stop(
        "\n  Got a score of NA for gene: ", curr_sym,
        "\n  expression_value: ", data_transform_compared$expression_value,
        "\n  pred_data_ref_expression: ", data_transform_compared$pred_data_ref_expression,
        "\n  with stretch_factor: ", stretch_factor,
        "\n  with curr_shift: ", curr_shift
      )
    }

    # TODO: should this be in attr() instead?
    all_scores[i] <- score
    all_data_transform_mean[i] <- data_transform_mean
    all_data_ref_mean[i] <- data_ref_mean
    all_data_transform_sd[i] <- data_transform_sd
    all_data_ref_sd[i] <- data_ref_sd
  }

  out <- data.table::data.table(
    data.frame(
      "gene" = curr_sym,
      "stretch" = stretch_factor,
      "shift" = all_shifts,
      "score" = all_scores,
      "data_transform_compared_mean" = all_data_transform_mean,
      "data_ref_compared_mean" = all_data_ref_mean,
      "data_transform_compared_sd" = all_data_transform_sd,
      "data_ref_compared_sd" = all_data_ref_sd
    )
  )

  return(out)
}
