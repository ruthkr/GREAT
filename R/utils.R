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
