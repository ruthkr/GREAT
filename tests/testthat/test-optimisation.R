brapa_sample_data <- data.table::fread(system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "greatR"))
reference <- "Ro18"
query <- "Col0"
gene_data <- brapa_sample_data[gene_id == "BRAA03G051930.3C"]
processed_data <- preprocess_data(gene_data, reference, query)
all_data <- processed_data$all_data

test_that("get_search_space_limits works", {
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

test_that("get_search_space_limits_from_params works", {
  stretches <- c(1, 2, 3)
  shifts <- c(-4, 4)
  space_lims <- get_search_space_limits_from_params(stretches = stretches, shifts = shifts)

  # Expected outputs
  expect_equal(names(space_lims), c("stretch_init", "stretch_lower", "stretch_upper", "shift_init", "shift_lower", "shift_upper"))
  expect_equal(space_lims$stretch_init, mean(stretches), tolerance = 1e-2)
  expect_equal(space_lims$stretch_lower, min(stretches), tolerance = 1e-2)
  expect_equal(space_lims$stretch_upper, max(stretches), tolerance = 1e-2)
  expect_equal(space_lims$shift_init, 0, tolerance = 1e-2)
  expect_equal(space_lims$shift_lower, min(shifts), tolerance = 1e-2)
  expect_equal(space_lims$shift_upper, max(shifts), tolerance = 1e-2)
})

test_that("calc_overlapping_percent works", {
  all_data_reg <- apply_registration(all_data, 2.75, 3.6)
  overlapping_raw <- calc_overlapping_percent(all_data)
  overlapping_reg <- calc_overlapping_percent(all_data_reg)

  # Expected outputs
  expect_gte(overlapping_reg, overlapping_raw)
  expect_equal(overlapping_reg, 1, tolerance = 1e-1)
})

test_that("objective_fun works", {
  # Expected outputs
  expect_equal(objective_fun(all_data, 2.75, 3.6, 0.5), -11.19, tolerance = 1e-2)
  expect_equal(objective_fun(all_data), -999)
})

test_that("optimise works", {
  optimisation_config = list(num_iterations = 1, num_fun_evals = 1)
  results_sa <- optimise(all_data, optimisation_config = optimisation_config, optimise_fun = optimise_using_sa)
  results_nm <- optimise(all_data, optimisation_config = optimisation_config, optimise_fun = optimise_using_nm)

  # Expected outputs
  expect_equal(names(results_sa), c("stretch", "shift", "loglik_score"))
  expect_equal(names(results_nm), c("stretch", "shift", "loglik_score"))
})
