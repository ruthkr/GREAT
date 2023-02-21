#' Compare curve dynamics for Hypothesis H2
#'
#' H2: the dasets are best explained by two different models.
#'
#' @noRd
calc_loglik_H2 <- function(data) {
  # Cut down to the data for each model
  data_query <- data[data$accession == "query"]
  data_ref <- data[data$accession == "ref"]

  # Fit using cubic spline with K+3 params
  num_spline_params <- 6

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

  return(loglik_separate)
}

#' Compare curve dynamics for Hypothesis H1
#'
#' H1: the dasets are best explained by one common model.
#'
#' @noRd
calc_loglik_H1 <- function(data) {
  # Cut down to the data for each model
  data_query <- data[data$accession == "query"]
  data_ref <- data[data$accession == "ref"]

  # Fit using cubic spline with K+3 params
  num_spline_params <- 6

  # Fit model for combined data together
  fit_all <- stats::lm(
    expression_value ~ splines::bs(timepoint, df = num_spline_params, degree = 3),
    data = data
  )

  # Calculate the log likelihoods and BIC
  loglik_ref <- calc_loglik(fit_all, data_ref)
  loglik_query <- calc_loglik(fit_all, data_query)
  loglik_combined <- loglik_query + loglik_ref

  return(loglik_combined)
}

#' Compare H1 and H2 hypotheses for registration
#'
#' @noRd
compare_H1_and_H2 <- function(data, stretch, shift, loglik_H1, loglik_H2) {
  # Specify parameters used to compute BIC
  num_spline_params <- 6
  num_registration_params <- 2
  num_obs <- nrow(data)

  # Calculate BIC for the models
  BIC_H1 <- calc_BIC(loglik_H1, num_spline_params + num_registration_params, num_obs)
  BIC_H2 <- calc_BIC(loglik_H2, 2 * num_spline_params, num_obs)

  model_comparison <- data.table::data.table(
    gene_id = unique(data$gene_id),
    BIC_separate = BIC_H2,
    BIC_combined = BIC_H1,
    stretch = stretch,
    shift = shift,
    registered = BIC_H2 > BIC_H1
  )

  return(model_comparison)
}
