#' Apply registration
#'
#' @param data Input data frame, either containing all replicates of gene expression or not.
#' @param stretches Candidate registration stretch factors to apply to query data.
#' @param shifts Candidate registration shift values to apply to query data.
#'
#' @noRd
apply_registration <- function(data, stretch, shift) {
  # Suppress "no visible binding for global variable" note
  accession <- NULL
  timepoint <- NULL
  time_delta <- NULL

  data <- data.table::copy(data)

  # Apply stretch
  data[, time_delta := as.numeric(time_delta)]
  data[, time_delta := if (accession == "query") time_delta * stretch else time_delta, by = .(accession)]

  # Apply shift
  data[, timepoint := as.numeric(timepoint)]
  data[, timepoint := if (accession == "query") time_delta + shift + min(timepoint) else timepoint, by = .(accession)]
  data[, c("time_delta") := NULL]

  return(data)
}
