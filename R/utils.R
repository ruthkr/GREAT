.onLoad <- function(libname, pkgname) {
  options(cli.progress_show_after = 0)
  options(cli.progress_clear = FALSE)
}

#' Validate names
#'
#' @noRd
match_names <- function(x, lookup) {
  unmatched <- x[-grep(paste(lookup, collapse = "$|"), x)]
  if (length(unmatched) > 0) {
    stop("Valid names are ", paste(lookup, collapse = ", "))
  }
}
