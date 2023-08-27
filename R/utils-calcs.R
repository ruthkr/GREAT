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
  # Read expression and variance
  expression_value <- data$expression_value
  sigma_squared <- data$var

  # Predict expressions
  pred_expression_value <- stats::predict(model, newdata = data)
  dist_squared <- (pred_expression_value - expression_value)^2

  # Calculate logLik
  loglik <- -sum(dist_squared / (2 * sigma_squared))

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

#' Calculate variance for the observed expression data
#'
#' @noRd
calc_variance <- function(all_data, exp_sd = NA) {
  # Suppress "no visible binding for global variable" note
  gene_id <- NULL
  accession <- NULL
  expression_value <- NULL
  var <- NULL

  if (any(is.na(exp_sd))) {
    # TODO: notify user
    # Calculate expression variance for replicates
    # all_data[, var := sd(expression_value)^2, by = .(gene_id, accession, timepoint)]

    # Calculate global expression variance
    all_data[, var := (diff(range(expression_value)) / 10)^2, by = .(gene_id, accession)]
  } else {
    all_data[, var := exp_sd^2]
  }

  return(all_data)
}
