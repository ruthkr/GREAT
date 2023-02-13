# # Get data where have expression data
# # and fit with a model
# sample_data <- system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "greatR") |>
#   utils::read.csv() |>
#   dplyr::select(
#     timepoint,
#     expression_value
#   ) |>
#   head(10)
#
# fit_model <- lm(timepoint ~ expression_value, data = sample_data)
#
# test_that("calc_loglik works", {
#   expect_equal(class(calc_loglik(fit_model, sample_data)), "numeric")
# })
#
#
# test_that("calc_BIC works", {
#   expect_equal(class(calc_BIC(calc_loglik(fit_model, sample_data), 2, 10)), "numeric")
# })
#
#
