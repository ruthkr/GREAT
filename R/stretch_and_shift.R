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

  return(data.table::data.table(data))
}
