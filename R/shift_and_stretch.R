#' Get the best result
#'
#' @noRd
get_best_result <- function(df) {
  # return TRUE/FALSE vector. TRUE for the smallest score
  # if tied for this, true for the one with the smallest stretch. (1x is smaller than 0.75x though)
  # if tied, then the one with the smallest shift
  is_best <- df$score == min(df$score)

  if (sum(is_best) == 1) {
    return(is_best)
  } else {
    cand_stretches <- df$stretch[is_best]
    # get the stretch with the best score, with the smallest divergence from 1
    min_stretch <- unique(cand_stretches[abs(cand_stretches - 1) == min(abs(cand_stretches - 1))])
    is_best[df$stretch != min_stretch] <- FALSE

    if (sum(is_best) == 1) {
      return(is_best)
    } else {
      cand_shifts <- df$shift[is_best]
      min_shift <- unique(cand_shifts[cand_shifts == min(cand_shifts)])
      is_best[df$shift != min_shift] <- FALSE
      if (sum(is_best) == 1) {
        return(is_best)
      } else {
        stop("error in get_best_result, somehow STILL more than one best shift tied?")
      }
    }
  }
}

#' Wrapper of applying best shifts and compare the registered and unregistered models
#'
#' @noRd
calculate_all_model_comparison_stats <- function(all_data_df,
                                                 best_shifts,
                                                 accession_data_to_transform,
                                                 accession_data_ref,
                                                 time_to_add) {
  if (!(accession_data_to_transform %in% unique(all_data_df$accession) & accession_data_ref %in% unique(all_data_df$accession))) {
    stop("error in calculate_all_model_comparison_stats() :
         all_data_df doesn't have the correct accession info - should have been
         converted to Ro18 & Col0")
  }

  # Apply the registration to the all rep data, so can use for model comparison
  shifted_all_data_df <- apply_best_shift(
    data = all_data_df,
    best_shifts,
    accession_data_to_transform,
    accession_data_ref,
    time_to_add
  )

  genes <- unique(shifted_all_data_df$locus_name)

  out.sepBIC <- numeric(length = length(genes))
  out.combBIC <- numeric(length = length(genes))

  i <- 0
  cli::cli_progress_step("Calculating registration vs non-registration comparison BIC ({i}/{length(genes)})", spinner = TRUE)
  for (i in seq_along(genes)) {
    curr_sym <- genes[i]

    BIC_comparison_list <- compare_registered_to_unregistered_model(
      curr_sym,
      original_data = all_data_df,
      data_df = shifted_all_data_df,
      accession_data_to_transform,
      accession_data_ref
    )

    out.sepBIC[i] <- BIC_comparison_list$separate.BIC
    out.combBIC[i] <- BIC_comparison_list$combined.BIC
    cli::cli_progress_update(force = TRUE)
  }

  model_comparison_dt <- data.table::data.table(
    "gene" = genes,
    "separate.BIC" = out.sepBIC,
    "registered.BIC" = out.combBIC
  )

  return(model_comparison_dt)
}

#' Register all expression over time using optimal shift found
#'
#' @param data Input data frame containing all replicates of gene expression in each genotype at each time point.
#' @param best_shifts Data frame of best shift parameters for each gene.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_ref Accession name of reference data.
#' @param time_to_add Time to be added when applying shift.
#'
#' @noRd
apply_best_shift <- function(data,
                             best_shifts,
                             accession_data_to_transform,
                             accession_data_ref,
                             time_to_add) {
  processed_data <- data.table::copy(data)

  processed_data <- apply_stretch(
    data = processed_data,
    best_shifts,
    accession_data_to_transform,
    accession_data_ref,
    time_to_add
  )

  # Normalise the expression data (If was normalised when calculating the expression data, is recorded in the _compared_mean, and _compared_sd columns. If no normalisation was carried out, then these should have values of 0 and 1. This was done using get_best_shift()).
  transform_was_not_normalised <- all(unique(best_shifts$data_transform_compared_mean) == 0)
  ref_was_not_normalised <- all(unique(best_shifts$data_ref_compared_mean) == 0)
  if (!transform_was_not_normalised | !ref_was_not_normalised) {
    processed_data <- apply_best_normalisation(
      data = processed_data,
      best_shifts,
      accession_data_to_transform,
      accession_data_ref
    )
  } else {
    # If no scaling carried out DURING the registration step
    cli::cli_alert_warning("No normalisation was carried out DURING registration (though may have been, prior to the comparison)")
  }

  # For each gene, shift the data to transform expression by the optimal shift found previously
  i <- 0
  cli::cli_progress_step("Applying best shift ({i}/{length(unique(processed_data$locus_name))})", spinner = TRUE)
  for (curr_gene in unique(processed_data$locus_name)) {
    curr_best_shift <- best_shifts$shift[best_shifts$gene == curr_gene]
    processed_data$shifted_time[processed_data$accession == accession_data_to_transform & processed_data$locus_name == curr_gene] <- processed_data$shifted_time[processed_data$accession == accession_data_to_transform & processed_data$locus_name == curr_gene] + curr_best_shift

    cli::cli_progress_update(force = TRUE)
    i <- i + 1
  }

  return(processed_data)
}

#' Apply stretch factor
#'
#' @noRd
apply_stretch <- function(data,
                          best_shifts,
                          accession_data_to_transform,
                          accession_data_ref,
                          time_to_add) {
  # Suppress "no visible binding for global variable" note
  delta_time <- NULL
  timepoint <- NULL
  accession <- NULL

  # Copy data
  data <- data.table::copy(data)

  # Stretch the expression of data to transform, leave reference data as is
  data[, delta_time := timepoint - min(timepoint), by = .(accession)]

  # Filter data based on the accession
  data_ref <- data[data$accession == accession_data_ref, ]
  data_to_transform <- data[data$accession == accession_data_to_transform, ]

  # Get the info of the strecth factor and merge data into one single data frame
  data_to_transform <- merge(
    data_to_transform,
    best_shifts[, c("gene", "stretch")],
    by.x = "locus_name",
    by.y = "gene"
  )

  data_to_transform$delta_time <- data_to_transform$delta_time * data_to_transform$stretch
  data_to_transform$stretch <- NULL

  # Bind by rows reference data and data to transform which have been stretched
  data <- rbind(data_ref, data_to_transform)

  # Record the stretched times (before individual shifting applied)
  # record the time (from start of timecourse) after stretching
  data$stretched_time_delta <- data$delta_time
  data$shifted_time <- data$delta_time

  # After stretching, add the time to the first datapoint back on
  data$shifted_time[data$accession == accession_data_to_transform] <- data$shifted_time[data$accession == accession_data_to_transform] + time_to_add
  data$shifted_time[data$accession == accession_data_ref] <- data$shifted_time[data$accession == accession_data_ref] + time_to_add
  data$delta_time <- NULL

  return(data)
}

#' Apply normalisation (after applying stretch)
#'
#' @noRd
apply_best_normalisation <- function(data,
                                     best_shifts,
                                     accession_data_to_transform,
                                     accession_data_ref) {
  i <- 0
  cli::cli_progress_step("Normalising expression by mean and sd of compared values ({i}/{length(unique(data$locus_name))})", spinner = TRUE)
  for (curr_gene in unique(data$locus_name)) {
    data_transform_mean <- best_shifts$data_transform_compared_mean[best_shifts$gene == curr_gene]
    data_ref_mean <- best_shifts$data_ref_compared_mean[best_shifts$gene == curr_gene]
    data_transform_sd <- best_shifts$data_transform_compared_sd[best_shifts$gene == curr_gene]
    data_ref_sd <- best_shifts$data_ref_compared_sd[best_shifts$gene == curr_gene]

    # If was compared
    if (length(data_transform_mean) != 0) {

      # Make sure that sd is not 0, since we do not want to divide by 0
      if (data_transform_sd != 0) {
        data$expression_value[data$locus_name == curr_gene & data$accession == accession_data_to_transform] <- (data$expression_value[data$locus_name == curr_gene & data$accession == accession_data_to_transform] - data_transform_mean) / data_transform_sd
      } else {
        data$expression_value[data$locus_name == curr_gene & data$accession == accession_data_to_transform] <- (data$expression_value[data$locus_name == curr_gene & data$accession == accession_data_to_transform] - data_transform_mean)
      }

      # Make sure that sd is not 0, since we do not want to divide by 0
      if (data_ref_sd != 0) {
        data$expression_value[data$locus_name == curr_gene & data$accession == accession_data_ref] <- (data$expression_value[data$locus_name == curr_gene & data$accession == accession_data_ref] - data_ref_mean) / data_ref_sd
      } else {
        data$expression_value[data$locus_name == curr_gene & data$accession == accession_data_ref] <- (data$expression_value[data$locus_name == curr_gene & data$accession == accession_data_ref] - data_ref_mean)
      }

      if (any(is.na(data$expression_value))) {
        stop("Have NAs in expression_value after rescaling in apply best_normalisation() for gene: {unique(data$locus_name)}")
      }
    } else {
      data$expression_value[data$locus_name == curr_gene & data$accession == accession_data_to_transform] <- NA
      data$expression_value[data$locus_name == curr_gene & data$accession == accession_data_ref] <- NA
    }

    cli::cli_progress_update(force = TRUE)
    i <- i + 1
  }

  return(data)
}

#' Comparing registered to unregistered model
#'
#' @noRd
compare_registered_to_unregistered_model <- function(curr_sym,
                                                     original_data,
                                                     data_df,
                                                     accession_data_to_transform,
                                                     accession_data_ref) {
  curr_data_df <- data_df[data_df$locus_name == curr_sym]

  # Flag the time points to be used in the modelling, only the ones which overlap!
  curr_data_df <- get_compared_timepoints(
    curr_data_df,
    accession_data_to_transform,
    accession_data_ref
  )

  # Cut down to the data for each model
  data_to_transform_spline <- curr_data_df[curr_data_df$is_compared == TRUE & curr_data_df$accession == accession_data_to_transform, ]
  data_ref_spline <- curr_data_df[curr_data_df$is_compared == TRUE & curr_data_df$accession == accession_data_ref, ]
  combined_spline_data <- curr_data_df[curr_data_df$is_compared == TRUE, ]


  # Fit the models - fit regression splines.
  # http://www.utstat.utoronto.ca/reid/sta450/feb23.pdf
  # For cubic spline, K+3 params where K=num.knots as can omit constant term
  num.spline.params <- 6 # number of parameters for each spline fitting (degree and this used to calculate num knots).
  num.registration.params <- 2 # stretch, shift
  num.obs <- nrow(combined_spline_data)

  data_to_transform_fit <- stats::lm(
    expression_value ~ splines::bs(shifted_time, df = num.spline.params, degree = 3),
    data = data_to_transform_spline
  )
  data_ref_fit <- stats::lm(
    expression_value ~ splines::bs(shifted_time, df = num.spline.params, degree = 3),
    data = data_ref_spline
  )
  combined_fit <- stats::lm(
    expression_value ~ splines::bs(shifted_time, df = num.spline.params, degree = 3),
    data = combined_spline_data
  )

  # Additional fit before transformation - using timepoint before registration
  # Get compared timepoint first:
  # ref_time <- curr_data_df %>%
  #   dplyr::filter(is_compared == TRUE, accession == accession_data_ref) %>%
  #   dplyr::pull(timepoint)
  #
  # trans_time <- curr_data_df %>%
  #   dplyr::filter(is_compared == TRUE, accession == accession_data_to_transform) %>%
  #   dplyr::pull(timepoint)

  # data_ref_spline_original <- original_data %>%
  #   dplyr::filter(timepoint %in% ref_time,
  #                 accession == accession_data_ref)
  #
  # data_to_transform_spline_original <- original_data %>%
  #   dplyr::filter(timepoint %in% trans_time,
  #                 accession == accession_data_to_transform)

  # data_ref_spline_original <- original_data %>%
  #   dplyr::filter(accession == accession_data_ref)
  #
  # data_to_transform_spline_original <- original_data %>%
  #   dplyr::filter(accession == accession_data_to_transform)
  #
  # data_ref_fit_before_reg <- stats::lm(
  #   expression_value ~ splines::bs(timepoint, df = num.spline.params, degree = 3),
  #   data = data_ref_spline_original
  # )
  #
  # data_to_transform_before_reg <- stats::lm(
  #   expression_value ~ splines::bs(timepoint, df = num.spline.params, degree = 3),
  #   data = data_to_transform_spline_original
  # )

  # Calculate the log likelihoods
  data_to_transform_logLik <- stats::logLik(data_to_transform_fit)
  data_ref_logLik <- stats::logLik(data_ref_fit)
  separate_logLik <- data_to_transform_logLik + data_ref_logLik # logLikelihoods, so sum
  combined_logLik <- stats::logLik(combined_fit)

  # Additional fit before transformation
  # data_to_transform_logLik_before_reg <- stats::logLik(data_to_transform_before_reg)
  # data_ref_logLik_before_reg <- stats::logLik(data_ref_fit_before_reg)
  # separate_logLik_before_reg <- data_to_transform_logLik_before_reg + data_ref_logLik_before_reg


  # Calculate the comparison.stats BIC: smaller is better!
  # 2*num.spline.params as fitting separate models for Ara * Col
  separate.BIC <- calc_BIC(separate_logLik, 2 * num.spline.params, num.obs)
  combined.BIC <- calc_BIC(combined_logLik, num.spline.params + num.registration.params, num.obs)
  # separate.BIC.before <- calc_BIC(separate_logLik_before_reg, 2 * num.spline.params, nrow(original_data))

  # Results object
  results_list <- list(
    separate.BIC = separate.BIC,
    combined.BIC = combined.BIC
  )

  return(results_list)
}
