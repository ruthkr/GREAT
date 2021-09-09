#' Registering data
#'
#' `scale_and_register_data` is a function to process all data. This includes scaling data before registration, finding and calculate score of optimal shifts and stretches, apply the best shifts and stretches.
#'
#' @param mean_df Input data frame contains the mean gene expression of each gene in each genotype at each timepoint.
#' @param all_data_df Input data frame contains all replicates of gene expression in each genotype at each timepoint.
#' @param stretches Candidate registration stretch factors to apply to data to transform.
#' @param shift_extreme The absolute maximum value which can be applied as a shift to gene expression timecourse (days).
#' @param num_shifts Number of shifts between minimum and maximum values of shift.
#' @param min_num_overlapping_points Number of minimum overlapping time points.  Shifts will be only considered if it leaves at least these many overlapping points after applying the registration function.
#' @param initial_rescale Scaling gene expression prior to registration if TRUE.
#' @param do_rescale Scaling gene expression using only overlapping timepoints points during registration.
#' @param testing Showing immediate results (including plots) if TRUE.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_fix Accession name of data fix.
#' @param data_to_transform_time_added Time points to be added in data to transform.
#' @param data_fix_time_added Time points to be added in data fix.
#'
#' @return List of dataframes: (a) `mean_df` is unchanged by `scale_and_register_data()`, (b) `mean_df_sc` is identical to `mean_df`, with additional column "sc.mean.cpm", (c) `imputed_mean_df` is registered expression data, (d) `all_shifts` is a table of candidate registraions applied, and score for each, and (e) `model_comparison_dt` is a table comparing the optimal registration function for each gene (based on `all_shifts` scores) to model with no registration applied.
#' @export
scale_and_register_data <- function(mean_df,
                                    all_data_df,
                                    stretches,
                                    shift_extreme,
                                    num_shifts,
                                    min_num_overlapping_points,
                                    initial_rescale,
                                    do_rescale,
                                    testing,
                                    accession_data_to_transform,
                                    accession_data_fix,
                                    data_to_transform_time_added,
                                    data_fix_time_added) {


  # Apply normalisation of expression for each gene across all timepoints
  mean_df_sc <- data.table::copy(mean_df)

  # Apply scaling
  mean_df_sc[, sc.mean_cpm := scale(mean_cpm, scale = TRUE, center = TRUE),
    by = .(locus_name, accession)]

  # Apply scaling before registration (if initial_rescale == TRUE), otherwise using original data
  if (initial_rescale == TRUE) {

    # apply rescale to mean_df prior to registration
    to_shift_df <- data.table::copy(mean_df_sc)
    to_shift_df$mean_cpm <- to_shift_df$sc.mean_cpm
    to_shift_df$sc.mean_cpm <- NULL

    # apply THE SAME rescale to all_data_df prior to registration
    all_data_df <- scale_all_rep_data(mean_df, all_data_df, "scale")
  } else {
    to_shift_df <- data.table::copy(mean_df)
  }

  message("Max value of mean_cpm of all_data_df :", max(all_data_df$mean_cpm))


  # calculate the best registration. Returns all tried registrations, best stretch and shift combo,
  # and AIC/BIC stats for comparison of best registration model to separate models for expression of
  # each gene in Ro18 and Col0.
  L <- get_best_stretch_and_shift(
    to_shift_df,
    all_data_df,
    stretches,
    do_rescale,
    min_num_overlapping_points,
    shift_extreme,
    num_shifts,
    testing,
    accession_data_to_transform,
    accession_data_fix,
    data_to_transform_time_added,
    data_fix_time_added
  )

  all_shifts <- L[["all_shifts"]]
  best_shifts <- L[["best_shifts"]]
  model_comparison_dt <- L[["model_comparison_dt"]]

  # browser()

  # message("Max value of all_shifts mean_cpm :", max(all_shifts$mean_cpm))

  # Add columns which flags which BIC and AIC values are better
  model_comparison_dt$BIC_registered_is_better <- (model_comparison_dt$registered.BIC < model_comparison_dt$separate.BIC)
  model_comparison_dt$AIC_registered_is_better <- (model_comparison_dt$registered.AIC < model_comparison_dt$separate.AIC)
  model_comparison_dt$ABIC_registered_is_better <- (model_comparison_dt$BIC_registered_is_better & model_comparison_dt$AIC_registered_is_better)


  # Report model comparison results
  print("################## Model comparison results #######################")
  print(paste0("AIC finds registration better than separate for :", sum(model_comparison_dt$AIC_registered_is_better), " / ", nrow(model_comparison_dt)))
  print(paste0("BIC finds registration better than separate for :", sum(model_comparison_dt$BIC_registered_is_better), " / ", nrow(model_comparison_dt)))
  print(paste0("AIC & BIC finds registration better than separate for :", sum(model_comparison_dt$ABIC_registered_is_better), " / ", nrow(model_comparison_dt)))
  print("###################################################################")


  # get the best-shifted and stretched mean gene expression, only to genes which registration is better than
  # separate models by BIC. Don't stretch out, or shift genes for which separate is better.
  # registration is applied to col0.
  shifted_mean_df <- apply_shift_to_registered_genes_only(to_shift_df,
                                                          best_shifts,
                                                          model_comparison_dt,
                                                          accession_data_to_transform,
                                                          accession_data_fix,
                                                          data_to_transform_time_added,
                                                          data_fix_time_added)

  message("Max value of mean_cpm :", max(shifted_mean_df$mean_cpm))

  # Impute transformed values at times == to the observed data fix points for each shifted transformed gene so can compare using heat maps.
  # transformed curves are the ones that been shifted around. Linear impute values for these
  # curves so that data fix samples can be compared to an transformed data point.
  imputed_mean_df <- impute_transformed_exp_values(shifted_mean_df,
                                               accession_data_to_transform,
                                               accession_data_fix)

  out <- list(
    "mean_df" = mean_df,
    "mean_df_sc" = mean_df_sc,
    "imputed_mean_df" = imputed_mean_df,
    "all_shifts" = all_shifts,
    "model_comparison_dt" = model_comparison_dt
  )

}


#' Scaling all un-averaged data
#'
#' `scale_all_rep_data` is a function to scale
#'
#' @param mean_df Input data containing mean of each time point.
#' @param all_rep_data Input all data (without taking mean).
#' @param scale_func Scaling method choice applied in all_rep_data. There are two options: (a) "scale" where all expression values are subtracted by mean value and divided by standard deviation and (b) "my_scale" where expression values are divided by mean values
#'
#' @return Scaled expression data in all_rep_data.
#' @export
scale_all_rep_data <- function(mean_df,
                               all_rep_data,
                               scale_func) {

  # apply the same scaling which done to the mean expression data
  # to all the reps.
  # (subtract mean, and divide by sd), using the values for the mean data
  # as this is what was used to find the best shift.

  # Calculate the summary stats to use for the rescaling
  gene_expression_stats <- unique(mean_df[, .(
    mean_val = mean(mean_cpm),
    sd_val = stats::sd(mean_cpm)
  ),
  by = .(locus_name, accession)
  ])

  # Combine all_rep_data with gene_expression_stats
  all_rep_data <- merge(all_rep_data, gene_expression_stats, by = c("locus_name", "accession"))

  # Adjust scaling calculation depends on the scale function choice
  if (scale_func == "scale") {
    all_rep_data$scaled_norm_cpm <- (all_rep_data$mean_cpm - all_rep_data$mean_val) / all_rep_data$sd_val
  } else if (scale_func == "my_scale") {
    all_rep_data$scaled_norm_cpm <- (all_rep_data$mean_cpm / all_rep_data$mean_val)
  } else {
    message("invalid scale option for scale_all_rep_data")

    stop()
  }

  out <- subset(all_rep_data,
    select = c(
      "locus_name",
      "accession",
      "tissue",
      "timepoint",
      "scaled_norm_cpm"
    )
  )

  names(out)[names(out) == "scaled_norm_cpm"] <- "mean_cpm"

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
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_fix Accession name of data fix.
#' @param data_to_transform_time_added Time points to be added in data to transform.
#' @param data_fix_time_added Time points to be added in data fix.
#'
#' @return List of data frames (a) all_shifts : all the combos of stretching and shifting tried for each gene, (b) best_shifts : the best stretch and shift combo found for each gene, as well as info for scaling, and (c) model_comparison.dt : AIC / BIC scores for best registerd model found, compared to separate model for each genes expression in the 2 accessions.
#' @export
get_best_stretch_and_shift <- function(to_shift_df,
                                       all_data_df,
                                       stretches,
                                       do_rescale,
                                       min_num_overlapping_points,
                                       shift_extreme,
                                       num_shifts,
                                       testing,
                                       accession_data_to_transform,
                                       accession_data_fix,
                                       data_to_transform_time_added,
                                       data_fix_time_added) {


  # Warning to make sure users have correct accession data
  if (!(accession_data_to_transform %in% all_data_df$accession & accession_data_fix %in% all_data_df$accession)) {
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
      testing,
      accession_data_to_transform,
      accession_data_fix
    )

    all_shifts <- unique(all_shifts) # ensure no duplicated rows

    # Cut down to single best shift for each gene
    all_shifts[, is_best:=get_best_result(.SD), by=.(gene)]
    best_shifts <- all_shifts[is_best==TRUE,]
    all_shifts$is_best <- NULL

    if (nrow(best_shifts)!= length(unique(all_data_df$locus_name))) {
      stop('get_best_stretch_and_shift() : got non-unique best shifts in best_shifts')
    }

    # Calculate the BIC & AIC for the best shifts found with this stretch.compared to treating the
    # gene's expression separately in data to transform and data fix
    model_comparison_dt <- calculate_all_model_comparison_stats(
      all_data_df,
      best_shifts,
      is_testing = testing,
      accession_data_to_transform,
      accession_data_fix,
      data_to_transform_time_added,
      data_fix_time_added
    )


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
  all_model_comparison_dt <- do.call('rbind', all_model_comparison_dt) # model comparison of best shift (for each stretch) to separate modeles

  # Get the best registration applied (best stretch, and best shift) for each gene,
  # picking by BIC alone will favour fewer overlapping (considered) data points.
  # Pick best in order to maximise how much better register.BIC is than separate.BIC
  all_model_comparison_dt$delta.BIC <- all_model_comparison_dt$registered.BIC - all_model_comparison_dt$separate.BIC

  # Best is one for which registered.BIC is as small as possible compared to separate.BIC
  all_model_comparison_dt[, is_best := (delta.BIC == min(delta.BIC)), by = .(gene)]
  best_model_comparison.dt <- all_model_comparison_dt[all_model_comparison_dt$is_best==TRUE]

  # If there is a tie for best registration for a gene, keep the first one as the best
  if (any(duplicated(best_model_comparison.dt$gene))) {
    print(paste0('found ', sum(duplicated(best_model_comparison.dt$gene)), ' tied optimal registrations. Removing duplicates'))
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


#' Apply shift for all registered genes
#'
#' `apply_shift_to_registered_genes_only` is a function to apply shift for all registered model based on `model_comparison_dt` using information from `best_shifts`.
#'
#' @param to_shift_df Input data frame.
#' @param best_shifts Data frame containing information of best shift and stretch values.
#' @param model_comparison_dt Data frame containing information of comparison of BIC and AIC for registred and non-registered genes.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_fix Accession name of data fix.
#' @param data_to_transform_time_added Time points to be added in data to transform.
#' @param data_fix_time_added Time points to be added in data fix.
#'
#' @return Data frame for all transformed genes for those with better BIC values.
#' @export
apply_shift_to_registered_genes_only <- function(to_shift_df,
                                                 best_shifts,
                                                 model_comparison_dt,
                                                 accession_data_to_transform,
                                                 accession_data_fix,
                                                 data_to_transform_time_added = 11,
                                                 data_fix_time_added) {

  # Genes for which registration model is better than separate model
  gene_to_register <- model_comparison_dt$gene[model_comparison_dt$BIC_registered_is_better]


  # Apply the registration transformation to these genes --------------------
  if (length(gene_to_register > 0)) {
    register.dt <- to_shift_df[to_shift_df$locus_name %in% gene_to_register, ]
    registered_dt <- apply_best_shift(
      data = register.dt,
      best_shifts,
      accession_data_to_transform,
      accession_data_fix,
      data_to_transform_time_added,
      data_fix_time_added
    )

    registered_dt$is_registered <- TRUE
  }

  # Genes for which the separate model is better than registration model
  genes_to_keep_separate <- model_comparison_dt$gene[!(model_comparison_dt$BIC_registered_is_better)]

  # Generate the columns for these needed to concat with registered_dt
  separate_dt <- to_shift_df[to_shift_df$locus_name %in% genes_to_keep_separate, ]
  # In order to ensure that separate copy
  separate_dt$stretched_time_delta <- 0

  # Apply the stretch transformation to these genes --------------------
  separate_dt[, stretched_time_delta := timepoint - min(timepoint), by = .(locus_name, accession)]

  # Here, we need to add additional time to make it comparable between data to transform and data fix
  # Therefore need to to this here, to keep unregistered in same frame as stretch 1, shift 0 registered genes.
  separate_dt$shifted_time <- separate_dt$stretched_time_delta + data_to_transform_time_added

  separate_dt$is_registered <- FALSE

  # Combine both registered and non-registered data frame
  if (length(gene_to_register > 0)) {
    out_dt <- rbind(registered_dt, separate_dt)
  } else {
    out_dt <- separate_dt
  }

  return(out_dt)
}



#' Setting transformed expression data and data fix to be the same in a set of common time points
#'
#' `impute_transformed_exp_values` is a function to impute transformed times at set of common time points in order to allow sample distance comparison to data fix. this means that transformed expression data were imputed relative to data fix time points. Since the original value of transformed data are not meant to be discarded, the imputed times are generated from minimum and maximum shifted time points of transformed data (not just data fix time points).
#'
#' @param shifted_mean_df All registered data frame.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_fix Accession name of data fix.
#'
#' @return
#' @export
impute_transformed_exp_values <- function(shifted_mean_df,
                                          accession_data_to_transform,
                                          accession_data_fix) {

  # The imputed transformed data times going to estimate gene expression for
  imputed_timepoints <- round(seq(min(shifted_mean_df$shifted_time), max(shifted_mean_df$shifted_time)))

  out_list <- list()
  out_list <- c(out_list, list(shifted_mean_df[shifted_mean_df$accession == accession_data_fix]))


  count <- 0
  for (curr_gene in unique(shifted_mean_df$locus_name)) {
    if (count %% 100 == 0) {
      message(count, " / ", length(unique(shifted_mean_df$locus_name)))
    }

    # Get the current gene expression data
    curr_df <- shifted_mean_df[shifted_mean_df$locus_name == curr_gene, ]

    transformed_df <- curr_df[curr_df$accession == accession_data_to_transform, ]
    # bra.df <- curr_df[curr_df$accession == accession_data_fix, ]

    interp_transformed_df <- data.table::data.table(
      "locus_name" = curr_gene,
      "accession" = accession_data_to_transform,
      "tissue" = "apex",
      "timepoint" = NA,
      "stretched_time_delta" = NA,
      "shifted_time" = imputed_timepoints,
      "is_registered" = unique(transformed_df$is_registered)[1]
    )

    # For each data fix timepoint, interpolate the comparable transformed expression data
    # by linear interpolation between the neighbouring two transformed expression values.
    # If not between two transformed expression values because shifted outside comparable range, set to NA.
    interp_transformed_df$mean_cpm <- sapply(imputed_timepoints,
      interpolate_data_fix_comparison_expression,
      data_fix_dt = transformed_df
    )

    out_list <- c(out_list, list(interp_transformed_df))
    count <- count + 1
  }

  out_df <- do.call("rbind", out_list)

  return(out_df)
}



#' @export
fix.accessions <- function(df, original.transformed.accession, original.other.accession) {

  # swap Col0 with original.transformed.accession, and Ro18 with original.other.accession
  new.df.accession <- df$accession
  new.df.accession[df$accession=='Col0'] <- original.transformed.accession
  new.df.accession[df$accession=='Ro18'] <- original.other.accession
  df$accession <- new.df.accession

  return(df)
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


