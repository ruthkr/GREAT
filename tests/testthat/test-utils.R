brapa_sample_data <- data.table::fread(system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "greatR"))

test_that("get_approximate_stretch works", {
  approximate_stretch <- get_approximate_stretch(
    brapa_sample_data,
    reference = "Ro18",
    query = "Col0"
  )

  # Expected outputs
  expect_equal(approximate_stretch, 2.6666666, tolerance = 1e-6)
  expect_gte(approximate_stretch, 0)
  expect_equal(class(approximate_stretch), "numeric")
})

test_that("match_names works", {
  a <- LETTERS[1:3]
  b <- LETTERS[4:5]

  # Expected outputs
  expect_error(match_names(x = a, lookup = b))
  expect_no_error(match_names(x = a, lookup = a))
})
