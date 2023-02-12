test_that("calc_BIC works", {
  expect_equal(class(calc_BIC(-3, 4, 1)), "numeric")
})
