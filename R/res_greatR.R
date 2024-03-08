new_res_greatR <- function(x) {
  structure(x, class = c("res_greatR", class(x)))
}

#' @export
print.res_greatR <- function(x, ...) {
  print(x$model_comparison)
  return(invisible(x))
}
