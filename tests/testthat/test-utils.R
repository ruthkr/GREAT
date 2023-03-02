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

test_that("validate_params works", {
  # Expected outputs
  expect_no_error(suppressMessages(validate_params(stretches = 1, shifts = 0, registration_type = "optimisation")))
  expect_no_error(suppressMessages(validate_params(stretches = NA, shifts = NA, registration_type = "optimisation")))
  expect_error(suppressMessages(validate_params(stretches = 1, shifts = NA, registration_type = "optimisation")))
  expect_no_error(suppressMessages(validate_params(stretches = 1, shifts = 0, registration_type = "manual")))
  expect_error(suppressMessages(validate_params(stretches = NA, shifts = NA, registration_type = "manual")))
})

test_that("cross_join works", {
  dt_a <- data.table(x = 1:2, y = 4:5)
  dt_b <- data.table(z = 1:3)
  dt_cj <- cross_join(dt_a, dt_b)

  # Expected outputs
  expect_equal(dim(dt_cj)[1], nrow(dt_a) * nrow(dt_b))
  expect_equal(dim(dt_cj)[2], ncol(dt_a) + ncol(dt_b))
  expect_equal(colnames(dt_cj), c(colnames(dt_a), colnames(dt_b)))
})
