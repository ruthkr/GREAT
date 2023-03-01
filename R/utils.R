.onLoad <- function(libname, pkgname) {
  options(cli.progress_show_after = 0)
  options(cli.progress_clear = FALSE)
}

#' Validate names
#'
#' @noRd
match_names <- function(x, lookup, error = NULL, name_string = "names", lookup_vec_last = " and ") {
  unmatched <- setdiff(x, lookup)
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
      stop(
        cli::format_error(c(
          "{.var stretches} and {.var shifts} must be {.cls numeric} vectors.",
          "x" = "You supplied vectors with {.cls NA} values."
        )),
        call. = FALSE
      )
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
