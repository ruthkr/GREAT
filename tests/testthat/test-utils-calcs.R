# Get expression data and fit with a model
sample_data <- data.table::fread(system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "greatR"))
sample_data_raw <- sample_data[1:10]
sample_data <- calc_variance(data.table::copy(sample_data_raw))

fit_model <- stats::lm(timepoint ~ expression_value, data = sample_data)

test_that("calc_loglik works", {
  loglik <- calc_loglik(fit_model, sample_data)

  # Expected outputs
  expect_equal(class(loglik), "numeric")
  expect_lte(loglik, 0)
  expect_equal(loglik, -163.775658, tolerance = 1e-6)
})

test_that("calc_BIC works", {
  loglik <- calc_loglik(fit_model, sample_data)
  bic <- calc_BIC(loglik, 2, 10)

  # Expected outputs
  expect_equal(class(bic), "numeric")
  expect_gte(bic, 0)
  expect_equal(bic, 332.156487, tolerance = 1e-6)
  expect_true(is.infinite(calc_BIC(loglik, 2, 0)))
  expect_warning(calc_BIC(loglik, 2, -10))
})

test_that("calc_variance works", {
  sample_data_with_reps <- data.table::copy(sample_data_raw)
  sample_data_with_no_reps <- data.table::copy(sample_data_raw[grep("-a", sample_data_raw$replicate), ])
  var_with_reps <- calc_variance(sample_data_with_reps)
  var_with_no_reps <- calc_variance(sample_data_with_no_reps)

  # Expected outputs
  expect_gt(length(unique(var_with_reps$var)), 1)
  expect_equal(length(unique(var_with_no_reps$var)), 1)
})
