brapa_sample_data <- data.table::fread(system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "greatR"))
reference <- "Ro18"
query <- "Col0"
gene_data <- brapa_sample_data[gene_id == "BRAA03G051930.3C"]
all_data <- preprocess_data(gene_data, reference, query, scaling_method = "z-score")

test_that("match_names works", {
  a <- LETTERS[1:3]
  b <- LETTERS[4:5]

  # Expected outputs
  expect_error(match_names(x = a, lookup = b))
  expect_error(match_names(x = a[1:2], lookup = a))
  expect_no_error(match_names(x = a, lookup = a))
})

test_that("validate_params works", {
  # Expected outputs
  expect_no_error(suppressMessages(validate_params(stretches = 1, shifts = 0, registration_type = "optimisation")))
  expect_no_error(suppressMessages(validate_params(stretches = NA, shifts = NA, registration_type = "optimisation")))
  expect_no_error(suppressMessages(validate_params(stretches = 1, shifts = NA, registration_type = "optimisation")))
  expect_no_error(suppressMessages(validate_params(stretches = 1, shifts = 0, registration_type = "manual")))
  expect_error(suppressMessages(validate_params(stretches = NA, shifts = NA, registration_type = "manual")))
  expect_error(suppressMessages(validate_params(stretches = 1, shifts = NA, registration_type = "manual")))
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

test_that("get_search_space_limits (auto) works", {
  space_lims <- get_search_space_limits(all_data)

  # Expected outputs
  expect_equal(names(space_lims), c("stretch_init", "stretch_lower", "stretch_upper", "shift_init", "shift_lower", "shift_upper"))
  expect_equal(space_lims$stretch_init, 2.667, tolerance = 1e-2)
  expect_equal(space_lims$stretch_lower, 1.333, tolerance = 1e-2)
  expect_equal(space_lims$stretch_upper, 4, tolerance = 1e-2)
  expect_equal(space_lims$shift_init, 0, tolerance = 1e-2)
  expect_equal(space_lims$shift_lower, -20, tolerance = 1e-2)
  expect_equal(space_lims$shift_upper, 16, tolerance = 1e-2)
  expect_no_error(get_search_space_limits(all_data))
  expect_error(get_search_space_limits(gene_data))
})

test_that("get_search_space_limits (bound) works", {
  stretches <- c(1, 2, 3)
  shifts <- c(-4, 4)
  space_lims <- get_search_space_limits(all_data, stretches = stretches, shifts = shifts)

  # Expected outputs
  expect_equal(names(space_lims), c("stretch_init", "stretch_lower", "stretch_upper", "shift_init", "shift_lower", "shift_upper"))
  expect_equal(space_lims$stretch_init, mean(stretches), tolerance = 1e-2)
  expect_equal(space_lims$stretch_lower, min(stretches), tolerance = 1e-2)
  expect_equal(space_lims$stretch_upper, max(stretches), tolerance = 1e-2)
  expect_equal(space_lims$shift_init, 0, tolerance = 1e-2)
  expect_equal(space_lims$shift_lower, min(shifts), tolerance = 1e-2)
  expect_equal(space_lims$shift_upper, max(shifts), tolerance = 1e-2)
})

test_that("get_search_space_limits (init) works", {
  stretches <- 1
  shifts <- 4
  space_lims <- get_search_space_limits(all_data, stretches = stretches, shifts = shifts)

  # Expected outputs
  expect_equal(names(space_lims), c("stretch_init", "stretch_lower", "stretch_upper", "shift_init", "shift_lower", "shift_upper"))
  expect_equal(space_lims$stretch_init, stretches, tolerance = 1e-2)
  expect_equal(space_lims$stretch_lower, 1.333, tolerance = 1e-2)
  expect_equal(space_lims$stretch_upper, 4, tolerance = 1e-2)
  expect_equal(space_lims$shift_init, shifts, tolerance = 1e-2)
  expect_equal(space_lims$shift_lower, -20, tolerance = 1e-2)
  expect_equal(space_lims$shift_upper, 16, tolerance = 1e-2)
})

test_that("calc_overlapping_percent works", {
  all_data_reg <- apply_registration(all_data, 3.10, 2.13)
  overlapping_raw <- calc_overlapping_percent(all_data)
  overlapping_reg <- calc_overlapping_percent(all_data_reg)

  # Expected outputs
  expect_gte(overlapping_reg, overlapping_raw)
  expect_equal(overlapping_reg, 1, tolerance = 1e-1)
})
