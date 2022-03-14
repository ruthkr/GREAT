test_that("calc_BIC works", {
  expect_equal(class(calc_BIC(-3, 4, 1)), "numeric")
})

test_that("calc_score works", {
  set.seed(1)
  data_to_transform_expression <- rnorm(5)
  data_ref_expression <- rnorm(5)
  score <- calc_score(data_to_transform_expression, data_ref_expression)
  expect_equal(class(score), "numeric")
  expect_gte(score, 0)
})
