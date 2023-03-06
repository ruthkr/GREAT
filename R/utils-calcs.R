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

#' Fit using cubic spline with K+3 parameters
#'
#' @param data Input data
#' @param x Predictor variable, by default "timepoint".
#' @param num_spline_params Number of parameters, or degrees of freedom, for each spline fitting. This is used to calculate the number of \code{knots}.
#' @param degree Degree of the piecewise polynomial, default is 3 for cubic splines.
#'
#' @noRd
fit_spline_model <- function(data, x = "timepoint", num_spline_params = 4, degree = 3) {
  fit_object <- stats::lm(
    stats::as.formula(
      paste("expression_value ~ splines::bs(", x, ", df = num_spline_params, degree = degree)")
    ),
    data = data
  )

  return(fit_object)
}
