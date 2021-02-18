#' Send message when running function
#'
#' @param fun Function name.
#'
#' @noRd
message_function_header <- function(fun) {
  # fun_message <- paste0(
  #   paste0(rep(" ", 3), collapse = ""),
  #   "● ",
  #   "RUNNING ",
  #   fun, "() ",
  #   "FUNCTION...",
  #   paste0(rep(" ", 3), collapse = "")
  # )
  # message(rep("-", times = nchar(fun_message)))
  # message(fun_message)
  # message(rep("-", times = nchar(fun_message)), "\n")
}

