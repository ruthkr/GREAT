#' Registering data
#'
#' `scale_and_register_data` is a function to process all data. This includes scaling data before registration, finding and calculate score of optimal shifts and stretches, apply the best shifts and stretches.
#'
#' @param input_df Input data frame contains all replicates of gene expression in each genotype at each timepoint.
#' @param stretches Candidate registration stretch factors to apply to data to transform.
#' @param shift_extreme The absolute maximum value which can be applied as a shift to gene expression timecourse (days).
#' @param num_shifts Number of shifts between minimum and maximum values of shift.
#' @param min_num_overlapping_points Number of minimum overlapping time points.  Shifts will be only considered if it leaves at least these many overlapping points after applying the registration function.
#' @param initial_rescale Scaling gene expression prior to registration if \code{TRUE}.
#' @param do_rescale Scaling gene expression using only overlapping timepoints points during registration.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_ref Accession name of reference data.
#' @param start_timepoint Start timepoint used to ... Time points to be added in both reference data and data to transform after shifting and stretching. Can be either \code{"reference"} (the default), \code{"transform"}, or \code{"zero"}.
#' @param expression_value_threshold Expression value threshold. Remove expressions if maximum is less than the threshold. If \code{NULL} keep all data.
#' @param is_data_normalised TRUE if dataset has been normalised prior to registration process.
#'
#' @return List of dataframes: (a) `mean_df` is unchanged by `scale_and_register_data()`, (b) `mean_df_sc` is identical to `mean_df`, with additional column `sc.expression_value`, (c) `imputed_mean_df` is registered expression data, (d) `all_shifts_df` is a table of candidate registrations applied, and score for each, and (e) `model_comparison_dt` is a table comparing the optimal registration function for each gene (based on `all_shifts_df` scores) to model with no registration applied.
#'
#' @export
scale_and_register_data <- function(input_df,
                                    stretches,
                                    shift_extreme,
                                    num_shifts,
                                    min_num_overlapping_points,
                                    initial_rescale,
                                    do_rescale,
                                    accession_data_to_transform,
                                    accession_data_ref,
                                    start_timepoint = c("reference", "transform", "zero"),
                                    expression_value_threshold = 5,
                                    is_data_normalised = FALSE) {
  # Validate parameters
  start_timepoint <- match.arg(start_timepoint)

  # Make sure the data are data.tables
  all_data_df <- data.table::as.data.table(input_df)

  mean_df <- get_mean_data(
    exp = all_data_df,
    expression_value_threshold = expression_value_threshold,
    accession_data_to_transform = accession_data_to_transform,
    is_data_normalised = is_data_normalised
  )

  # Parse start_timepoint
  if (start_timepoint == "reference") {
    time_to_add <- min(all_data_df[accession == accession_data_ref, timepoint])
  } else if (start_timepoint == "transform") {
    time_to_add <- min(all_data_df[accession == accession_data_to_transform, timepoint])
  } else {
    time_to_add <- 0
  }

  # Filter genes of original input data as in mean_df
  all_data_df <- all_data_df[all_data_df$locus_name %in% unique(mean_df$locus_name)]
  all_data_df <- subset(all_data_df, select = c("locus_name", "accession", "tissue", "timepoint", "expression_value", "group"))

  # TODO: validate colnames
  # Apply normalisation of expression for each gene across all timepoints
  mean_df_sc <- data.table::copy(mean_df)

  # Apply scaling
  mean_df_sc[, sc.expression_value := scale(expression_value, scale = TRUE, center = TRUE), by = .(locus_name, accession)]

  # Apply scaling before registration (if initial_rescale == TRUE), otherwise using original data
  if (initial_rescale == TRUE) {

    # apply rescale to mean_df prior to registration
    to_shift_df <- data.table::copy(mean_df_sc)
    to_shift_df$expression_value <- to_shift_df$sc.expression_value
    to_shift_df$sc.expression_value <- NULL

    # apply THE SAME rescale to all_data_df prior to registration
    all_data_df <- scale_all_rep_data(mean_df, all_data_df, "scale")
  } else {
    to_shift_df <- data.table::copy(mean_df)
  }

  cli::cli_h1("Information before registration")
  cli::cli_alert_info("Max value of expression_value of all_data_df: {cli::col_cyan(round(max(all_data_df$expression_value), 2))}")

  # Calculate the best registration. Returns all tried registrations, best stretch and shift combo,and AIC/BIC stats for comparison of best registration model to separate models for expression ofeach gene in Ro18 and Col0
  cli::cli_h1("Analysing models for all stretch and shift factor")

  L <- get_best_stretch_and_shift(
    to_shift_df,
    all_data_df,
    stretches,
    do_rescale,
    min_num_overlapping_points,
    shift_extreme,
    num_shifts,
    accession_data_to_transform,
    accession_data_ref,
    time_to_add
  )

  all_shifts <- L[["all_shifts"]]
  best_shifts <- L[["best_shifts"]]
  model_comparison_dt <- L[["model_comparison_dt"]]

  # Add columns which flags which BIC and AIC values are better
  model_comparison_dt$BIC_registered_is_better <- (model_comparison_dt$registered.BIC < model_comparison_dt$separate.BIC)
  model_comparison_dt$AIC_registered_is_better <- (model_comparison_dt$registered.AIC < model_comparison_dt$separate.AIC)
  model_comparison_dt$ABIC_registered_is_better <- (model_comparison_dt$BIC_registered_is_better & model_comparison_dt$AIC_registered_is_better)


  # Report model comparison results
  cli::cli_h1("Model comparison results")
  cli::cli_alert_info("AIC finds registration better than separate for: {cli::col_cyan(sum(model_comparison_dt$AIC_registered_is_better), '/', nrow(model_comparison_dt))}")
  cli::cli_alert_info("BIC finds registration better than separate for: {cli::col_cyan(sum(model_comparison_dt$BIC_registered_is_better), '/', nrow(model_comparison_dt))}")
  cli::cli_alert_info("AIC & BIC finds registration better than separate for: {cli::col_cyan(sum(model_comparison_dt$ABIC_registered_is_better), '/', nrow(model_comparison_dt))}")

  # Get the best-shifted and stretched mean gene expression, only to genes which registration is better than
  # separate models by BIC. Don't stretch out, or shift genes for which separate is better.

  cli::cli_h1("Applying the best-shifts and stretches to gene expression")
  shifted_mean_df <- apply_shift_to_registered_genes_only(
    to_shift_df,
    best_shifts,
    model_comparison_dt,
    accession_data_to_transform,
    accession_data_ref,
    time_to_add
  )

  cli::cli_alert_info("Max value of expression_value: {cli::col_cyan(round(max(shifted_mean_df$expression_value), 2))}")

  # Impute transformed values at times == to the observed reference data points for each shifted transformed gene so can compare using heat maps.
  # Transformed curves are the ones that been shifted around. Linear impute values for these curves so that reference data samples can be compared to an transformed data point.
  imputed_mean_df <- impute_transformed_exp_values(
    shifted_mean_df,
    accession_data_to_transform,
    accession_data_ref
  )

  out <- list(
    "mean_df" = mean_df,
    "mean_df_sc" = mean_df_sc,
    "imputed_mean_df" = imputed_mean_df,
    "all_shifts_df" = all_shifts,
    "model_comparison_dt" = model_comparison_dt
  )
}

#' Scaling all un-averaged data
#'
#' `scale_all_rep_data` is a function to apply the same scaling which done to the mean expression data to all the reps. (subtract mean, and divide by sd), using the values for the mean data as this is what was used to find the best shift.
#'
#' @param mean_df Input data containing mean of each time point.
#' @param all_rep_data Input all data (without taking mean).
#' @param scale_func Scaling method choice applied in all_rep_data. There are two options: (a) "scale" where all expression values are subtracted by mean value and divided by standard deviation and (b) "my_scale" where expression values are divided by mean values.
#'
#' @return Scaled expression data in all_rep_data.
scale_all_rep_data <- function(mean_df,
                               all_rep_data,
                               scale_func) {
  # Calculate the summary statistics to use for the rescaling
  gene_expression_stats <- unique(
    mean_df[, .(
      mean_val = mean(expression_value),
      sd_val = stats::sd(expression_value)
    ), by = .(locus_name, accession)]
  )

  # Combine all_rep_data with gene_expression_stats
  all_rep_data <- merge(
    all_rep_data,
    gene_expression_stats,
    by = c("locus_name", "accession")
  )

  # Adjust scaling calculation depends on the scale function choice
  if (scale_func == "scale") {
    all_rep_data$scaled_norm_expression_value <- (all_rep_data$expression_value - all_rep_data$mean_val) / all_rep_data$sd_val
  } else if (scale_func == "my_scale") {
    all_rep_data$scaled_norm_expression_value <- (all_rep_data$expression_value / all_rep_data$mean_val)
  } else {
    stop("invalid scale option for scale_all_rep_data")
  }

  out <- subset(
    all_rep_data,
    select = c(
      "locus_name",
      "accession",
      "tissue",
      "timepoint",
      "scaled_norm_expression_value"
    )
  )

  names(out)[names(out) == "scaled_norm_expression_value"] <- "expression_value"

  return(out)
}

#' Calculate best shifts and stretches for each gene, also calculate AIC/BIC under registration or non-registration
#'
#' `get_best_stretch_and_shift` is a function to stretch in all stretches and calculates best shift, by comparing SUM of squares difference. For the best shift in each stretch, compares to separate models to calculate AIC/BIC under registration or no registration.
#'
#' @param to_shift_df Input data containing mean of each time point.
#' @param all_data_df Input all data (without taking mean).
#' @param stretches Vector data of stretches.
#' @param do_rescale Apply "scale" to compared points for each shift if \code{TRUE}, use original mean expression data if \code{FALSE}.
#' @param min_num_overlapping_points Bound the extreme allowed shifts, such than at least this many timepoints are being compared for both accessions.
#' @param shift_extreme Approximation of maximum and minimum shifts allowed.
#' @param num_shifts Number of different shifts to be considered.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_ref Accession name of reference data.
#' @param time_to_add Time points to be added in both reference data and data to transform after shifting and stretching.
#'
#' @return List of data frames (a) all_shifts : all the combos of stretching and shifting tried for each gene, (b) best_shifts : the best stretch and shift combo found for each gene, as well as info for scaling, and (c) model_comparison.dt : AIC / BIC scores for best registerd model found, compared to separate model for each genes expression in the 2 accessions.
#'
#' @export
get_best_stretch_and_shift <- function(to_shift_df,
                                       all_data_df,
                                       stretches,
                                       do_rescale,
                                       min_num_overlapping_points,
                                       shift_extreme,
                                       num_shifts,
                                       accession_data_to_transform,
                                       accession_data_ref,
                                       time_to_add) {
  # Warning to make sure users have correct accession data
  if (!(accession_data_to_transform %in% all_data_df$accession & accession_data_ref %in% all_data_df$accession)) {
    stop("get_best_stretch_and_shift(): data accessions should have been converted to correct accession.")
  }

  all_all_shifts <- rep(list(0), length(stretches))
  all_best_shifts <- rep(list(0), length(stretches))
  all_model_comparison_dt <- rep(list(0), length(stretches))

  for (i in 1:length(stretches)) {
    stretch <- stretches[i]
    cli::cli_h2("Analysing models for stretch factor = {stretch}")

    # Calculate all the shift scores given this stretch. Score is mean(dist^2), over overlapping points if do_rescale=T, is rescaled by the mean FOR THE OVERLAPPING POINTS. (but not by the SD.)
    all_shifts <- calculate_all_best_shifts(
      num_shifts,
      mean_df = to_shift_df,
      stretch_factor = stretch,
      do_rescale,
      shift_extreme,
      min_num_overlapping_points,
      accession_data_to_transform,
      accession_data_ref
    )

    all_shifts <- unique(all_shifts) # ensure no duplicated rows

    # Cut down to single best shift for each gene
    all_shifts[, is_best := get_best_result(.SD), by = .(gene)]
    best_shifts <- all_shifts[is_best == TRUE, ]
    all_shifts$is_best <- NULL

    if (nrow(best_shifts) != length(unique(all_data_df$locus_name))) {
      stop("get_best_stretch_and_shift(): got non-unique best shifts in best_shifts")
    }

    # Calculate the BIC & AIC for the best shifts found with this stretch.compared to treating the gene's expression separately in data to transform and reference data
    model_comparison_dt <- calculate_all_model_comparison_stats(
      all_data_df,
      best_shifts,
      accession_data_to_transform,
      accession_data_ref,
      time_to_add
    )

    # Add info on the stretch and shift applied
    model_comparison_dt <- merge(
      model_comparison_dt,
      best_shifts[, c("gene", "stretch", "shift"), ],
      by = "gene"
    )

    # Record the results for the current stretch factor
    all_all_shifts[[i]] <- all_shifts
    all_best_shifts[[i]] <- best_shifts
    all_model_comparison_dt[[i]] <- model_comparison_dt
    cli::cli_alert_success("Finished analysing models for stretch factor = {stretch}")
  }

  # all the combinations of shift, and stretch tried
  all_shifts <- do.call("rbind", all_all_shifts)
  # the best shifts for each stretch
  all_best_shifts <- do.call("rbind", all_best_shifts)
  # model comparison of best shift (for each stretch) to separate models
  all_model_comparison_dt <- do.call("rbind", all_model_comparison_dt)

  # Correct -Inf BIC and AIC values to -9999 so that delta.BIC is not Inf or NaN
  all_model_comparison_dt <- all_model_comparison_dt %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      dplyr::across(
        .cols = c(registered.BIC, registered.AIC, separate.BIC, separate.AIC),
        .fns = function(x) {
          if (!is.finite(x)) {
            x <- 9999 * sign(x)
          }
          return(x)
        }
      )
    ) %>%
    dplyr::ungroup() %>%
    data.table::as.data.table()

  # Get the best registration applied (best stretch, and best shift) for each gene, picking by BIC alone will favour fewer overlapping (considered) data points.
  # Pick best in order to maximise how much better register.BIC is than separate.BIC
  all_model_comparison_dt$delta.BIC <- all_model_comparison_dt$registered.BIC - all_model_comparison_dt$separate.BIC

  # Best is one for which registered.BIC is as small as possible compared to separate.BIC
  all_model_comparison_dt[, is_best := (delta.BIC == min(delta.BIC)), by = .(gene)]
  best_model_comparison.dt <- all_model_comparison_dt[all_model_comparison_dt$is_best == TRUE]

  # If there is a tie for best registration for a gene, keep the first one as the best
  if (any(duplicated(best_model_comparison.dt$gene))) {
    message("found ", sum(duplicated(best_model_comparison.dt$gene)), " tied optimal registrations. Removing duplicates")
    best_model_comparison.dt <- best_model_comparison.dt[!(duplicated(best_model_comparison.dt$gene)), ]
  }

  best_model_comparison.dt$delta.BIC <- NULL

  # Cut down best shifts to the best shift for the best stretch only
  best_shifts <- merge(
    all_best_shifts,
    best_model_comparison.dt[, c("gene", "stretch", "shift")],
    by = c("gene", "stretch", "shift")
  )

  # There should be only 1 best shift for each gene, stop if it is not the case
  if (!(nrow(best_shifts) == length(unique(to_shift_df$locus_name)))) {
    stop()
  }

  return(list(
    "all_shifts" = all_shifts,
    "best_shifts" = best_shifts,
    "model_comparison_dt" = best_model_comparison.dt
  ))
}

#' Apply shift for all registered genes
#'
#' `apply_shift_to_registered_genes_only` is a function to apply shift for all registered model based on `model_comparison_dt` using information from `best_shifts`.
#'
#' @param to_shift_df Input data frame.
#' @param best_shifts Data frame containing information of best shift and stretch values.
#' @param model_comparison_dt Data frame containing information of comparison of BIC and AIC for registred and non-registered genes.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_ref Accession name of reference data.
#' @param time_to_add Time points to be added in both reference data and data to transform after shifting and stretching.
#'
#' @return Data frame for all transformed genes for those with better BIC values.
apply_shift_to_registered_genes_only <- function(to_shift_df,
                                                 best_shifts,
                                                 model_comparison_dt,
                                                 accession_data_to_transform,
                                                 accession_data_ref,
                                                 time_to_add) {
  # Genes for which registration model is better than separate model
  gene_to_register <- model_comparison_dt$gene[model_comparison_dt$BIC_registered_is_better]

  # Apply the registration transformation to these genes
  if (length(gene_to_register > 0)) {
    register.dt <- to_shift_df[to_shift_df$locus_name %in% gene_to_register, ]
    registered_dt <- apply_best_shift(
      data = register.dt,
      best_shifts,
      accession_data_to_transform,
      accession_data_ref,
      time_to_add
    )

    registered_dt$is_registered <- TRUE
  }

  # Genes for which the separate model is better than registration model
  genes_to_keep_separate <- model_comparison_dt$gene[!(model_comparison_dt$BIC_registered_is_better)]

  if (length(genes_to_keep_separate) > 0) {
    # Generate the columns for these needed to concat with registered_dt
    separate_dt <- to_shift_df[to_shift_df$locus_name %in% genes_to_keep_separate, ]
    # In order to ensure that separate copy
    separate_dt$stretched_time_delta <- 0
    # Apply the stretch transformation to these genes
    separate_dt[, stretched_time_delta := timepoint - min(timepoint), by = .(locus_name, accession)]
    # Here, we need to add additional time to make it comparable between data to transform and reference data
    # Therefore need to to this here, to keep unregistered in same frame as stretch 1, shift 0 registered genes
    separate_dt$shifted_time <- separate_dt$stretched_time_delta + time_to_add
    separate_dt$is_registered <- FALSE
  }

  # Combine both registered and non-registered data frame
  if (length(gene_to_register) > 0 & length(genes_to_keep_separate) > 0) {
    out_dt <- rbind(registered_dt, separate_dt)
  } else if (length(gene_to_register) > 0 & length(genes_to_keep_separate) == 0) {
    out_dt <- registered_dt
  } else {
    out_dt <- separate_dt
  }

  return(out_dt)
}

#' Setting transformed expression data and reference data to be the same in a set of common time points
#'
#' `impute_transformed_exp_values` is a function to impute transformed times at set of common time points in order to allow sample distance comparison to reference data. this means that transformed expression data were imputed relative to reference data time points. Since the original value of transformed data are not meant to be discarded, the imputed times are generated from minimum and maximum shifted time points of transformed data (not just reference data time points).
#'
#' @param shifted_mean_df All registered data frame.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_ref Accession name of reference data.
impute_transformed_exp_values <- function(shifted_mean_df,
                                          accession_data_to_transform,
                                          accession_data_ref) {
  # The imputed transformed data times going to estimate gene expression for
  imputed_timepoints <- round(seq(min(shifted_mean_df$shifted_time), max(shifted_mean_df$shifted_time)))

  out_list <- list()
  out_list <- c(out_list, list(shifted_mean_df[shifted_mean_df$accession == accession_data_ref]))

  i <- 0
  cli::cli_progress_step("Imputing transformed expression values ({i}/{length(unique(shifted_mean_df$locus_name))})", spinner = TRUE)
  for (curr_gene in unique(shifted_mean_df$locus_name)) {
    # Get the current gene expression data
    curr_df <- shifted_mean_df[shifted_mean_df$locus_name == curr_gene, ]

    transformed_df <- curr_df[curr_df$accession == accession_data_to_transform, ]

    interp_transformed_df <- data.table::data.table(
      "locus_name" = curr_gene,
      "accession" = accession_data_to_transform,
      "tissue" = transformed_df$tissue[1],
      "timepoint" = NA,
      "stretched_time_delta" = NA,
      "shifted_time" = imputed_timepoints,
      "is_registered" = unique(transformed_df$is_registered)[1]
    )

    # For each reference data timepoint, interpolate the comparable transformed expression data by linear interpolation between the neighbouring two transformed expression values.
    # If not between two transformed expression values because shifted outside comparable range, set to NA.
    interp_transformed_df$expression_value <- sapply(
      imputed_timepoints,
      interpolate_data_ref_comparison_expression,
      data_ref_dt = transformed_df
    )

    out_list <- c(out_list, list(interp_transformed_df))

    cli::cli_progress_update(force = TRUE)
    i <- i + 1
  }

  out_df <- do.call("rbind", out_list)

  return(out_df)
}
