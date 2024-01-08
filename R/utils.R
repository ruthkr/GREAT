.onLoad <- function(libname, pkgname) {
  options(cli.progress_show_after = 0)
  options(cli.progress_clear = FALSE)
}

#' Validate names
#'
#' @noRd
match_names <- function(x, lookup, error = NULL, name_string = "names", lookup_vec_last = " and ") {
  unmatched <- c(setdiff(x, lookup), setdiff(lookup, x))
  if (length(unmatched) > 0) {
    stop(
      cli::format_error(c(
        error,
        "i" = "Valid {name_string} are {.val {cli::cli_vec(lookup, style = list(vec_sep = ', ', vec_last = {lookup_vec_last}))}}.",
        "x" = "You supplied {.val {x}}."
      )),
      call. = FALSE
    )
  }
}

#' Validate registration parameters
#'
#' @noRd
validate_params <- function(stretches, shifts, registration_type = c("optimisation", "manual")) {
  # Registration with optimisation
  if (registration_type == "optimisation") {
    if (all(is.na(stretches), is.na(shifts))) {
      cli::cli_alert_info("Using computed stretches and shifts search space limits.")
    } else if (all(!is.na(stretches), !is.na(shifts))) {
      cli::cli_alert_info("Using provided stretches and shifts to define search space limits.")
    } else {
      var_na <- ifelse(is.na(stretches), "stretches", "shifts")
      var_num <- ifelse(is.na(stretches), "shifts", "stretches")

      cli::cli_alert_info("Using provided {var_num} and computed {var_na} to define search space limits.")
    }
  }

  # Manual registration
  if (registration_type == "manual") {
    if (any(is.na(stretches), is.na(shifts))) {
      stop(
        cli::format_error(c(
          "{.var stretches} and {.var shifts} must be {.cls numeric} vectors.",
          "x" = "You supplied vectors with {.cls NA} values."
        )),
        call. = FALSE
      )
    }
  }
}

#' Perform crossing in {data.table}
#'
#' @noRd
cross_join <- function(a, b) {
  cj <- data.table::CJ(
    seq_len(nrow(a)),
    seq_len(nrow(b))
  )
  cbind(a[cj[[1]], ], b[cj[[2]], ])
}

#' Get approximate stretch factor
#'
#' \code{get_approximate_stretch()} is a function to get a stretch factor
#' estimation given input data. This function will take the time point ranges
#' of both reference and query data and compare them to estimate the stretch
#' factor.
#'
#' @param data Input data frame, either containing all replicates of gene expression or not.
#' @param reference Accession name of reference data.
#' @param query Accession name of query data.
#'
#' @return This function returns an estimation of a stretch factor for registering the data.
#'
#' @export
get_approximate_stretch <- function(data, reference = "ref", query = "query") {
  # Suppress "no visible binding for global variable" note
  accession <- NULL
  timepoint <- NULL
  time_range <- NULL

  # Make sure the data are data.tables
  data <- data.table::as.data.table(data)

  # Validate accession names
  match_names(
    x = c(reference, query),
    lookup = unique(data$accession),
    error = "Must review the supplied {.var reference} and {.var query} values:",
    name_string = "accession values"
  )

  # Calculate approximate stretch factor
  deltas <- data[, .(time_range = max(timepoint) - min(timepoint)), by = .(accession)]

  stretch_factor <- deltas[accession == reference, time_range] / deltas[accession == query, time_range]

  return(stretch_factor)
}

#' Calculate limits of the search space
#'
#' @noRd
get_search_space_limits <- function(data, stretches = NA, shifts = NA, overlapping_percent = 0.5) {
  stretch_space_lims <- get_stretch_search_space_limits(data, stretches)
  shift_space_lims <- get_shift_search_space_limits(data, shifts, stretch_space_lims$stretch_upper, overlapping_percent)
  space_lims <- c(stretch_space_lims, shift_space_lims)

  return(space_lims)
}

#' Calculate limits of the stretch search space
#'
#' @noRd
get_stretch_search_space_limits <- function(data, stretches = NA) {
  # Suppress "no visible binding for global variable" note
  accession <- NULL
  timepoint <- NULL

  # Check calculation mode
  if (all(is.na(stretches))) {
    calc_mode <- "auto"
  } else if (length(stretches) == 1) {
    calc_mode <- "init"
  } else if (length(stretches) >= 2) {
    calc_mode <- "bound"
  }

  # Calculate boundary
  if (calc_mode == "bound") {
    # Calculate limits
    stretch_init <- mean(stretches)
    stretch_lower <- min(stretches)
    stretch_upper <- max(stretches)
  } else {
    # Initial approximate value
    stretch_approx <- get_approximate_stretch(data)

    # Initial value
    if (calc_mode == "auto") {
      stretch_init <- stretch_approx
    } else if (calc_mode == "init") {
      stretch_init <- stretches
    }

    # Calculate limits
    stretch_lower <- 0.5 * stretch_approx
    stretch_upper <- 1.5 * stretch_approx
  }

  # Results object
  results_list <- list(
    stretch_init = stretch_init,
    stretch_lower = stretch_lower,
    stretch_upper = stretch_upper
  )

  return(results_list)
}

#' Calculate limits of the shift search space
#'
#' @noRd
get_shift_search_space_limits <- function(data, shifts = NA, stretch_upper, overlapping_percent = 0.5) {
  # Suppress "no visible binding for global variable" note
  accession <- NULL
  timepoint <- NULL

  # Check calculation mode
  if (all(is.na(shifts))) {
    calc_mode <- "auto"
  } else if (length(shifts) == 1) {
    calc_mode <- "init"
  } else if (length(shifts) >= 2) {
    calc_mode <- "bound"
  }

  # Calculate boundary
  if (calc_mode == "bound") {
    # Calculate limits
    shift_lower <- min(shifts)
    shift_upper <- max(shifts)
  } else {
    # Extract time point ranges
    timepoints_ref <- unique(data[accession == "ref", timepoint])
    timepoints_query <- unique(data[accession == "query", timepoint])

    # Calculate time point ranges
    range_ref <- diff(range(timepoints_ref))
    range_query <- diff(range(timepoints_query))
    range_query_max_stretch <- stretch_upper * range_query

    # Calculate minimum and maximum timepoints in which the curves overlap
    min_timepoint <- min(timepoints_ref) + overlapping_percent * range_ref - range_query_max_stretch
    max_timepoint <- max(timepoints_ref) - overlapping_percent * range_ref + range_query_max_stretch

    # Calculate limits
    shift_lower <- min_timepoint - min(timepoints_query)
    shift_upper <- (max_timepoint - range_query_max_stretch) - min(timepoints_query)
  }

  # Calculate initial value (zero if possible)
  if (calc_mode %in% c("auto", "bound")) {
    shift_init <- 0
  } else {
    shift_init <- shifts
  }

  if (shift_init < shift_lower | shift_init > shift_upper) {
    shift_init <- mean(c(shift_lower, shift_upper))
  }

  # Results object
  results_list <- list(
    shift_init = shift_init,
    shift_lower = shift_lower,
    shift_upper = shift_upper
  )

  return(results_list)
}

#' Calculate overlapping percentage between reference and query data time point ranges
#'
#' @noRd
calc_overlapping_percent <- function(data) {
  # Suppress "no visible binding for global variable" note
  accession <- NULL
  timepoint <- NULL

  # Extract time point ranges
  range_ref <- range(unique(data[accession == "ref", timepoint]))
  range_query <- range(unique(data[accession == "query", timepoint]))

  if (all(range_ref[2] >= range_query[2], range_ref[1] <= range_query[1])) {
    # Query is fully contained on reference
    overlapping_percent <- 1
  } else {
    # Calculate overlapping percent over reference
    overlap <- min(c(range_ref[2], range_query[2])) - max(c(range_ref[1], range_query[1]))
    overlapping_percent <- overlap / diff(range_ref)
  }

  return(overlapping_percent)
}
