#' Calculate Bayesian Information Criterion
#'
#' @noRd
calc_BIC <- function(logL, num_params, num_obs) {
  return((-2 * logL) + log(num_obs) * num_params)
}

#' Calculate log likelihood of the fitted model to the observed expression data
#'
#' @noRd
calc_loglik <- function(model, data) {
  y <- data$expression_value
  n <- length(y)
  preds <- stats::predict(model, newdata = data)
  dist_squared <- sum((preds - y)^2)
  sigma_squared <- sum((y - mean(y))^2, na.rm = TRUE) / (n - 1)

  # Calculate logLik
  loglik <- -dist_squared / (2 * sigma_squared)

  return(loglik)
}
