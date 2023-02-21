# Get expression data and fit with a model
sample_data <- data.table::fread(system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "greatR"))
sample_data <- sample_data[1:10, .(timepoint, expression_value)]

fit_model <- stats::lm(timepoint ~ expression_value, data = sample_data)

test_that("calc_loglik works", {
  loglik <- calc_loglik(fit_model, sample_data)

  # Expected outputs
  expect_equal(class(loglik), "numeric")
  expect_lte(loglik, 0)
  expect_equal(loglik, -3.777704, tolerance = 1e-6)
})

test_that("calc_BIC works", {
  loglik <- calc_loglik(fit_model, sample_data)
  bic <- calc_BIC(loglik, 2, 10)

  # Expected outputs
  expect_equal(class(bic), "numeric")
  expect_gte(bic, 0)
  expect_equal(bic, 12.160578, tolerance = 1e-6)
  expect_true(is.infinite(calc_BIC(loglik, 2, 0)))
  expect_warning(calc_BIC(loglik, 2, -10))
})
