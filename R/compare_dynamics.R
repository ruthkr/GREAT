#' Compare curve dynamics for Hypothesis H2
#'
#' H2: the dasets are best explained by two different models.
#'
#' @noRd
compare_dynamics_H2 <- function(data) {
  # Cut down to the data for each model
  data_query <- data[data$accession == "query"]
  data_ref <- data[data$accession == "ref"]

  # Fit using cubic spline with K+3 params where K=num.knots as can omit constant term
  num_spline_params <- 6 # number of parameters for each spline fitting (degree and this used to calculate num knots)
  num_obs <- nrow(data)

  # Fit models for reference and query
  fit_ref <- stats::lm(
    expression_value ~ splines::bs(timepoint, df = num_spline_params, degree = 3),
    data = data_ref
  )
  fit_query <- stats::lm(
    expression_value ~ splines::bs(timepoint, df = num_spline_params, degree = 3),
    data = data_query
  )

  # Calculate the log likelihoods and BIC
  loglik_query <- calc_loglik(fit_query, data_query)
  loglik_ref <- calc_loglik(fit_ref, data_ref)
  loglik_separate <- loglik_query + loglik_ref

  BIC <- calc_BIC(loglik_separate, 2 * num_spline_params, num_obs)

  return(BIC)
}

#' Compare curve dynamics for Hypothesis H1
#'
#' H1: the dasets are best explained by one common model.
#'
#' @noRd
compare_dynamics_H1 <- function(data) {
  # Cut down to the data for each model
  data_query <- data[data$accession == "query"]
  data_ref <- data[data$accession == "ref"]

  # Fit using cubic spline with K+3 params where K=num.knots as can omit constant term
  num_spline_params <- 6 # number of parameters for each spline fitting (degree and this used to calculate num knots)
  num_registration_params <- 2 # stretch, shift
  num_obs <- nrow(data)

  # Fit model for combined data together
  fit_all <- stats::lm(
    expression_value ~ splines::bs(timepoint, df = num_spline_params, degree = 3),
    data = data
  )

  # Calculate the log likelihoods and BIC
  loglik_ref <- calc_loglik(fit_all, data_ref)
  loglik_query <- calc_loglik(fit_all, data_query)
  loglik_combined <- loglik_query + loglik_ref

  BIC <- calc_BIC(loglik_combined, num_spline_params + num_registration_params, num_obs)

  return(BIC)
}

#' Compare H1 and H2 hypotheses for registration
#'
#' @noRd
compare_H1_and_H2 <- function(data, BIC_H1, BIC_H2) {
  model_comparison <- data.table::data.table(
    gene_id = unique(data$gene_id),
    BIC_separate = BIC_H2,
    BIC_combined = BIC_H1,
    stretch = attr(data, "stretch"),
    shift = attr(data, "shift"),
    registered = BIC_H2 > BIC_H1
  )

  return(model_comparison)
}
