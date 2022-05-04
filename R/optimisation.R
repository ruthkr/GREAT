#' Calculate best shifts and stretches for each gene, also calculate BIC under registration or non-registration
#'
#' Simplified version of \code{\link{get_best_stretch_and_shift}} for \code{\link{optimise_registration_params}}.
#'
#' @noRd
get_best_stretch_and_shift_simplified <- function(to_shift_df,
                                                  all_data_df,
                                                  stretches,
                                                  shifts,
                                                  do_rescale,
                                                  min_num_overlapping_points,
                                                  maintain_min_num_overlapping_points,
                                                  accession_data_to_transform,
                                                  accession_data_ref,
                                                  time_to_add) {
  # Warning to make sure users have correct accession data
  if (!(accession_data_to_transform %in% all_data_df$accession & accession_data_ref %in% all_data_df$accession)) {
    stop("get_best_stretch_and_shift(): data accessions should have been converted to correct accession.")
  }

  # Calculate all the shift scores given this stretch. Score is mean(dist^2), over overlapping points if do_rescale=T, is rescaled by the mean FOR THE OVERLAPPING POINTS. (but not by the SD.)
  all_shifts <- calculate_all_best_shifts(
    mean_df = to_shift_df,
    stretch_factor = stretches,
    shifts,
    do_rescale,
    min_num_overlapping_points,
    maintain_min_num_overlapping_points,
    accession_data_to_transform,
    accession_data_ref
  )

  # Ensure no duplicated rows
  all_shifts <- unique(all_shifts)

  # Calculate the BIC for the best shifts found with this stretch.compared to treating the gene's expression separately in data to transform and reference data
  model_comparison_dt <- calculate_all_model_comparison_stats(
    all_data_df,
    all_shifts,
    accession_data_to_transform,
    accession_data_ref,
    time_to_add
  )

  # Results object
  results_list <- list(
    model_comparison_dt = model_comparison_dt
  )

  return(results_list)
}


#' Get BIC score from registering data
#'
#' Simplified version of \code{\link{scale_and_register_data}} for \code{\link{optimise_registration_params}}.
#'
#' @noRd
get_BIC_from_registering_data <- function(input_df,
                                          stretches,
                                          shifts,
                                          min_num_overlapping_points,
                                          maintain_min_num_overlapping_points = FALSE,
                                          initial_rescale,
                                          do_rescale,
                                          accession_data_to_transform,
                                          accession_data_ref,
                                          start_timepoint = c("reference", "transform", "zero"),
                                          expression_value_threshold = 5,
                                          is_data_normalised = FALSE) {
  # Validate parameters
  start_timepoint <- match.arg(start_timepoint)

  # Preprocess data
  processed_data <- preprocess_data(
    input_df = input_df,
    initial_rescale = initial_rescale,
    accession_data_to_transform = accession_data_to_transform,
    accession_data_ref = accession_data_ref,
    start_timepoint = start_timepoint,
    expression_value_threshold = expression_value_threshold,
    is_data_normalised = is_data_normalised
  )

  all_data_df <- processed_data$all_data_df
  to_shift_df <- processed_data$to_shift_df
  time_to_add <- processed_data$time_to_add

  # Calculate the best registration
  best_registration_list <- get_best_stretch_and_shift_simplified(
    to_shift_df,
    all_data_df,
    stretches,
    shifts,
    do_rescale,
    min_num_overlapping_points,
    maintain_min_num_overlapping_points,
    accession_data_to_transform,
    accession_data_ref,
    time_to_add
  )

  registered_BIC <- best_registration_list$model_comparison_dt$registered.BIC
  separate_BIC <- best_registration_list$model_comparison_dt$separate.BIC
  BIC_diff <- registered_BIC - separate_BIC

  return(BIC_diff)
}

#' Calculate boundary box for Simulated Annealing
#'
#' @noRd
get_boundary_box <- function(input_df,
                             stretches_bound = NA,
                             shifts_bound = NA,
                             accession_data_to_transform,
                             accession_data_ref,
                             min_num_overlapping_points,
                             maintain_min_num_overlapping_points = FALSE,
                             expression_value_threshold) {
  # Validate parameters
  are_bounds_defined <- !any(c(any(is.na(stretches_bound)), any(is.na(shifts_bound))))

  # Initial stretch value
  stretch_init <- get_approximate_stretch(
    input_df = input_df,
    accession_data_to_transform = accession_data_to_transform,
    accession_data_ref = accession_data_ref
  )

  # User defined boundaries
  if (are_bounds_defined) {
    # Stretch limits
    cli::cli_alert_info("Using user-defined stretches and shifts as boundaries")

    if (length(stretches_bound) < 2) {
      bound_factor <- 0.5
      stretches_bound <- c(stretches_bound - bound_factor, stretches_bound + bound_factor)
    }

    stretch_lower <- min(stretches_bound)
    stretch_upper <- max(stretches_bound)

    if (length(shifts_bound) < 2) {
      bound_factor <- 0.5
      shifts_bound <- c(shifts_bound - bound_factor, shifts_bound + bound_factor)
    }

    shift_lower <- min(shifts_bound)
    shift_upper <- max(shifts_bound)

    shift_init <- mean(c(shift_lower, shift_upper))
  }

  # Computed boundaries
  if (!are_bounds_defined) {
    # Initial stretch limits
    cli::cli_alert_info("Using computed stretches and shifts boundaries")
    stretch_lower <- 0.5 * stretch_init
    stretch_upper <- 1.5 * ceiling(stretch_init)

    # Define stretch and shift limits given min_num_overlapping_points
    all_data_df <- data.table::as.data.table(input_df)

    mean_df <- get_mean_data(
      exp = all_data_df,
      expression_value_threshold = expression_value_threshold,
      accession_data_to_transform = accession_data_to_transform,
      is_data_normalised = FALSE
    )

    # Create candidate shifts limits data frame
    bound_limits <- seq(round(stretch_lower, 1), round(stretch_upper, 1), 0.1) %>%
      purrr::map(
        function(stretch) {
          extreme_shifts <- tryCatch(
            {
              shifts <- get_extreme_shifts_for_all(
                mean_df,
                stretch_factor = stretch,
                min_num_overlapping_points = min_num_overlapping_points,
                shift_extreme = 1000,
                accession_data_to_transform = accession_data_to_transform,
                accession_data_ref = accession_data_ref
              ) %>%
                suppressWarnings()

              df <- data.frame(
                stretch = stretch,
                shift_lower = shifts[[1]],
                shift_upper = shifts[[2]]
              )

              return(df)
            },
            error = function(error_message) {
              df <- data.frame(
                stretch = stretch,
                shift_lower = NA,
                shift_upper = NA
              )

              return(df)
            }
          )
        }
      ) %>%
      purrr::reduce(dplyr::bind_rows) %>%
      dplyr::filter(
        !is.na(shift_lower),
        !is.na(shift_upper),
        !is.infinite(shift_lower),
        !is.infinite(shift_upper)
      )

    if (!maintain_min_num_overlapping_points) {
      # Consider biggest box possible
      shift_upper <- max(bound_limits$shift_upper)
      shift_lower <- min(bound_limits$shift_lower)
    } else {
      # Filter bound_limits df only for those shift values between Q25% - Q75%
      quan <- unname(stats::quantile(bound_limits$shift_upper))

      bound_limits <- bound_limits %>%
        dplyr::filter(
          shift_upper <= quan[4],
          shift_upper >= quan[2]
        )

      # Restrict shift limits to make sure min_num_overlapping_points condition is maintained
      shift_upper <- min(bound_limits$shift_upper)
      shift_lower <- max(bound_limits$shift_lower)
    }

    # Define stretch limits
    stretch_lower <- min(bound_limits$stretch)
    stretch_upper <- max(bound_limits$stretch)
  }

  # Correct initial stretch value
  if (stretch_init < stretch_lower | stretch_init > stretch_upper) {
    stretch_init <- mean(c(stretch_lower, stretch_upper))
  }

  # Correct initial shift value
  shift_init <- 0
  if (shift_init < shift_lower | shift_init > shift_upper) {
    shift_init <- mean(c(shift_lower, shift_upper))
  }

  # Results object
  results_list <- list(
    stretch_init = stretch_init,
    stretch_lower = stretch_lower,
    stretch_upper = stretch_upper,
    shift_init = shift_init,
    shift_lower = shift_lower,
    shift_upper = shift_upper
  )

  return(results_list)
}

#' Optimise registration parameters with Simulated Annealing for single gene
#' @noRd
#' @importFrom rlang .data
optimise_registration_params_single_gene <- function(input_df,
                                                     initial_guess = NA,
                                                     stretches_bound = NA,
                                                     shifts_bound = NA,
                                                     initial_rescale = FALSE,
                                                     do_rescale = TRUE,
                                                     min_num_overlapping_points = 4,
                                                     maintain_min_num_overlapping_points = FALSE,
                                                     expression_value_threshold = 5,
                                                     accession_data_to_transform,
                                                     accession_data_ref,
                                                     start_timepoint,
                                                     is_data_normalised,
                                                     num_iterations = 60) {
  # Function to optimise
  BIC_diff <- function(x) {
    stretch <- x[1]
    shift <- x[2]

    tryCatch(
      {
        BIC <- get_BIC_from_registering_data(
          input_df = input_df,
          stretches = stretch,
          shifts = shift,
          min_num_overlapping_points = min_num_overlapping_points,
          maintain_min_num_overlapping_points = FALSE,
          expression_value_threshold = expression_value_threshold,
          initial_rescale = initial_rescale,
          do_rescale = do_rescale,
          accession_data_to_transform = accession_data_to_transform,
          accession_data_ref = accession_data_ref,
          start_timepoint = start_timepoint,
          is_data_normalised = is_data_normalised
        ) %>%
          suppressMessages() %>%
          suppressWarnings()

        return(BIC)
      },
      error = function(error_message) {
        BIC <- 999

        return(BIC)
      }
    )
  }

  # Calculate boundary box and initial guess
  boundary_box <- get_boundary_box(
    input_df,
    stretches_bound,
    shifts_bound,
    accession_data_to_transform,
    accession_data_ref,
    min_num_overlapping_points,
    maintain_min_num_overlapping_points,
    expression_value_threshold
  )

  if (any(is.na(initial_guess))) {
    stretch_init <- boundary_box$stretch_init
    shift_init <- boundary_box$shift_init
  } else {
    stretch_init <- initial_guess[1]
    shift_init <- initial_guess[2]
  }

  stretch_lower <- boundary_box$stretch_lower
  stretch_upper <- boundary_box$stretch_upper
  shift_lower <- boundary_box$shift_lower
  shift_upper <- boundary_box$shift_upper

  # Calculate cooling schedule
  t0 <- 1000
  t_min <- 0.1
  r_cooling <- (t_min / t0)^(1 / num_iterations)
  # TODO: Explore best default
  num_inner_loop_iter <- 100

  # Perform SA using {optimization}
  optim_sa_res <- optimization::optim_sa(
    fun = BIC_diff,
    start = c(stretch_init, shift_init),
    trace = TRUE,
    lower = c(stretch_lower, shift_lower),
    upper = c(stretch_upper, shift_upper),
    control = list(
      t0 = t0,
      t_min = t_min,
      nlimit = num_inner_loop_iter,
      r = r_cooling,
      rf = 3,
      dyn_rf = FALSE
    )
  )

  # Parse results
  locus_name <- unique(input_df$locus_name)

  result_df <- data.frame(
    gene = locus_name,
    stretch = round(optim_sa_res$par[1], 3),
    shift = round(optim_sa_res$par[2], 3),
    BIC_diff = optim_sa_res$function_value,
    is_registered = optim_sa_res$function_value < 0
  )

  trace_df <- optim_sa_res$trace %>%
    as.data.frame() %>%
    dplyr::mutate(
      gene = locus_name,
      is_registered = .data$loss < 0
    ) %>%
    dplyr::filter(.data$is_registered) %>%
    dplyr::select(
      .data$gene,
      stretch = .data$x_1,
      shift = .data$x_2,
      BIC_diff = .data$loss,
      .data$is_registered
    ) %>%
    dplyr::distinct()

  # Results object
  results_list <- list(
    optimum_params_df = result_df,
    candidate_params_df = trace_df
  )

  return(results_list)
}

#' Optimise registration parameters with Simulated Annealing
#'
#' @param input_df Input data frame containing all replicates of gene expression in each genotype at each time point.
#' @param genes List of genes to optimise.
#' @param stretches_bound Optional candidate registration stretch factors define search space, otherwise automatic.
#' @param shifts_bound Optional candidate registration shift values to define search space, otherwise automatic.
#' @param initial_rescale Scaling gene expression prior to registration if \code{TRUE}.
#' @param do_rescale Scaling gene expression using only overlapping time points points during registration.
#' @param min_num_overlapping_points Number of minimum overlapping time points. Shifts will be only considered if it leaves at least these many overlapping points after applying the registration function.
#' @param maintain_min_num_overlapping_points Whether to automatically calculate extreme (minimum and maximum) values of \code{shifts} to maintain specified \code{min_num_overlapping_points} condition. By default, \code{FALSE}.
#' @param accession_data_to_transform Accession name of data which will be transformed.
#' @param accession_data_ref Accession name of reference data.
#' @param start_timepoint Time points to be added in both reference data and data to transform after shifting and stretching. Can be either \code{"reference"} (the default), \code{"transform"}, or \code{"zero"}.
#' @param expression_value_threshold Expression value threshold. Remove expressions if maximum is less than the threshold. If \code{NULL} keep all data.
#' @param is_data_normalised \code{TRUE} if dataset has been normalised prior to registration process.
#' @param num_iterations Maximum number of iterations of the algorithm. Default is 100.
#'
#' @return List of optimum registration parameters, \code{optimum_params_df}, and other candidate registration parameters, \code{candidate_params_df} for all genes.
#' @export
#' @importFrom rlang .data
optimise_registration_params <- function(input_df,
                                         genes = NULL,
                                         stretches_bound = NA,
                                         shifts_bound = NA,
                                         initial_rescale = FALSE,
                                         do_rescale = TRUE,
                                         min_num_overlapping_points = 4,
                                         maintain_min_num_overlapping_points = FALSE,
                                         accession_data_to_transform,
                                         accession_data_ref,
                                         start_timepoint = c("reference", "transform", "zero"),
                                         expression_value_threshold = 5,
                                         is_data_normalised = FALSE,
                                         num_iterations = 60) {
  # Validate genes
  if (is.null(genes)) {
    genes <- unique(input_df$locus_name)
  }

  # Apply optimise_registration_params_single_gene() over all genes
  raw_results <- cli::cli_progress_along(
    genes,
    format = "{cli::pb_spin} Optimising registration parameters for genes ({cli::pb_current}/{cli::pb_total})",
    format_done = "{cli::col_green(cli::symbol$tick)} Optimising registration parameters for genes ({cli::pb_total}/{cli::pb_total}) {cli::col_white(paste0('[', cli::pb_elapsed, ']'))}",
    clear = FALSE
  ) %>%
    purrr::map(
      function(gene) {
        curr_df <- input_df %>%
          dplyr::filter(.data$locus_name == genes[[gene]])

        opt_res <- optimise_registration_params_single_gene(
          input_df = curr_df,
          stretches_bound = stretches_bound,
          shifts_bound = shifts_bound,
          initial_guess = NA,
          initial_rescale,
          do_rescale,
          min_num_overlapping_points,
          maintain_min_num_overlapping_points,
          expression_value_threshold,
          accession_data_to_transform,
          accession_data_ref,
          start_timepoint,
          is_data_normalised,
          num_iterations
        )

        return(opt_res)
      }
    )

  # Parse raw results
  optimum_params_reduced <- raw_results %>%
    purrr::map(~ purrr::pluck(.x, "optimum_params_df")) %>%
    purrr::reduce(dplyr::bind_rows)

  candidate_params_reduced <- raw_results %>%
    purrr::map(~ purrr::pluck(.x, "candidate_params_df")) %>%
    purrr::reduce(dplyr::bind_rows)

  # Results object
  results_list <- list(
    optimum_params_df = optimum_params_reduced,
    candidate_params_df = candidate_params_reduced
  )

  return(results_list)
}

#' Calculate best shifts and stretches for each gene, also calculate BIC under registration or non-registration
#'
#' @noRd
get_best_stretch_and_shift_after_optimisation <- function(to_shift_df,
                                                          all_data_df,
                                                          optimised_parameters,
                                                          do_rescale,
                                                          min_num_overlapping_points,
                                                          accession_data_to_transform,
                                                          accession_data_ref,
                                                          time_to_add) {
  # Suppress "no visible binding for global variable" note
  is_best <- NULL
  gene <- NULL
  delta.BIC <- NULL
  locus_name <- NULL

  # Warning to make sure users have correct accession data
  if (!(accession_data_to_transform %in% all_data_df$accession & accession_data_ref %in% all_data_df$accession)) {
    stop("get_best_stretch_and_shift(): data accessions should have been converted to correct accession.")
  }

  params_df <- optimised_parameters$optimum_params_df
  gene_list <- params_df$gene

  all_all_shifts <- rep(list(0), length(gene_list))
  all_best_shifts <- rep(list(0), length(gene_list))
  all_model_comparison_dt <- rep(list(0), length(gene_list))

  for (i in 1:length(gene_list)) {
    gene <- gene_list[i]
    gene_data_df <- all_data_df[locus_name == gene]
    cli::cli_h2("Analysing models for gene = {gene}")

    # Calculate all the shift scores given this stretch. Score is mean(dist^2), over overlapping points if do_rescale=T, is rescaled by the mean FOR THE OVERLAPPING POINTS. (but not by the SD.)
    all_shifts <- calculate_all_best_shifts(
      mean_df = to_shift_df[locus_name == gene],
      stretch_factor = params_df[params_df$gene == gene, ]$stretch,
      shifts = params_df[params_df$gene == gene, ]$shift,
      do_rescale,
      min_num_overlapping_points,
      maintain_min_num_overlapping_points = FALSE,
      accession_data_to_transform,
      accession_data_ref
    )

    # Ensure no duplicated rows
    all_shifts <- unique(all_shifts)

    # Cut down to single best shift for each gene (Alex's original logic)
    best_shifts <- all_shifts
    best_shifts$is_best <- TRUE

    if (nrow(best_shifts) != length(unique(gene_data_df$locus_name))) {
      stop("get_best_stretch_and_shift(): got non-unique best shifts in best_shifts")
    }

    # Calculate the BIC for the best shifts found with this stretch.compared to treating the gene's expression separately in data to transform and reference data
    model_comparison_dt <- calculate_all_model_comparison_stats(
      gene_data_df,
      best_shifts,
      accession_data_to_transform,
      accession_data_ref,
      time_to_add
    )

    # Add info on the stretch and shift applied
    model_comparison_dt <- merge(
      model_comparison_dt,
      best_shifts[, c("gene", "stretch", "shift", "score"), ],
      by = "gene"
    )

    # Record the results for the current stretch factor
    all_all_shifts[[i]] <- all_shifts
    all_best_shifts[[i]] <- best_shifts
    all_model_comparison_dt[[i]] <- model_comparison_dt
    cli::cli_alert_success("Finished analysing models for gene = {gene}")
  }

  # all the combinations of shift, and stretch tried
  all_shifts <- do.call("rbind", all_all_shifts)
  # the best shifts for each stretch
  all_best_shifts <- do.call("rbind", all_best_shifts)
  # model comparison of best shift (for each gene) to separate models
  all_model_comparison_dt <- do.call("rbind", all_model_comparison_dt)

  # Correct -Inf BIC values to -9999 so that delta.BIC is not Inf or NaN
  all_model_comparison_dt <- all_model_comparison_dt %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      dplyr::across(
        .cols = c(.data$registered.BIC, .data$separate.BIC),
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

  best_model_comparison.dt$BIC_registered_is_better <- best_model_comparison.dt$delta.BIC < 0
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

  # Results object
  results_list <- list(
    "all_shifts" = all_shifts,
    "best_shifts" = best_shifts,
    "model_comparison_dt" = best_model_comparison.dt
  )

  return(results_list)
}
