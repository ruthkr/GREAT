brapa_sample_data <- data.table::fread(system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "greatR"))
reference <- "Ro18"
query <- "Col0"
gene_data <- brapa_sample_data[gene_id == "BRAA03G051930.3C"]

# Preprocessing and intermediate functions ----

test_that("preprocess_data works", {
  processed_data <- preprocess_data(brapa_sample_data, reference, query, scaling_method = "scale")
  processed_data_norm <- preprocess_data(brapa_sample_data, reference, query, scaling_method = "normalise")
  all_data <- processed_data$all_data
  all_data_norm <- processed_data_norm$all_data

  # Expected outputs
  expect_equal(names(processed_data), "all_data")
  expect_equal(names(processed_data_norm), "all_data")
  expect_equal(class(all_data)[1], "data.table")
  expect_equal(class(all_data_norm)[1], "data.table")
  expect_gte(mean(all_data$expression_value), mean(all_data_norm$expression_value))
  expect_equal(colnames(all_data), c(colnames(brapa_sample_data), "time_delta"))
  expect_equal(levels(unique(all_data$accession)), c("ref", "query"))
  expect_equal(nrow(all_data), nrow(brapa_sample_data))
})

test_that("register_manually works", {
  gene_data <- preprocess_data(gene_data, reference, query)$all_data
  stretch <- 2.75
  shift <- 3.6
  loglik_separate <- -10.404
  results <- register_manually(gene_data, stretch, shift, loglik_separate)
  results_simple <- register_manually(gene_data, stretch, shift, loglik_separate, return_data_reg = FALSE)

  # Expected outputs
  expect_equal(names(results), c("data_reg", "model_comparison"))
  expect_equal(names(results_simple), "model_comparison")
  expect_equal(colnames(results$data_reg), c("gene_id", "accession", "timepoint", "replicate", "expression_value"))
  expect_equal(colnames(results$model_comparison), c("gene_id", "stretch", "shift", "BIC_separate", "BIC_combined", "registered"))
  expect_equal(results$model_comparison$stretch, stretch)
  expect_equal(results$model_comparison$shift, shift)
  expect_equal(results$model_comparison$registered, TRUE)
  expect_equal(results$model_comparison, results_simple$model_comparison)
})

# Full pipeline ----

test_that("register (with no optimisation) works", {
  stretch <- 2.75
  shift <- 3.6
  registration_results <- register(
    gene_data,
    reference = "Ro18",
    query = "Col0",
    stretches = stretch,
    shifts = shift,
    optimise_registration_parameters = FALSE
  ) |>
    suppressMessages()

  data_reg <- registration_results$data
  model_comparison <- registration_results$model_comparison

  # Expected outputs
  expect_equal(names(registration_results), c("data", "model_comparison"))
  expect_equal(colnames(data_reg), c("gene_id", "accession", "expression_value", "replicate", "timepoint", "timepoint_reg"))
  expect_equal(colnames(model_comparison), c("gene_id", "stretch", "shift", "BIC_separate", "BIC_combined", "registered"))
  expect_equal(model_comparison$registered, TRUE)
  expect_error(suppressMessages(register(gene_data, reference = "Ro18", query = "Col0", stretches = stretch, shifts = NA, optimise_registration_parameters = FALSE)))
})

test_that("register (with optimisation) works", {
  registration_results <- register(
    gene_data,
    reference = "Ro18",
    query = "Col0",
    optimisation_method = "nm",
    optimisation_config = list(num_iterations = 10, num_fun_evals = 10)
  ) |>
    suppressMessages()

  data_reg <- registration_results$data
  model_comparison <- registration_results$model_comparison

  # Expected outputs
  expect_equal(names(registration_results), c("data", "model_comparison"))
  expect_equal(colnames(data_reg), c("gene_id", "accession", "expression_value", "replicate", "timepoint", "timepoint_reg"))
  expect_equal(colnames(model_comparison), c("gene_id", "stretch", "shift", "BIC_separate", "BIC_combined", "registered"))
  expect_equal(model_comparison$registered, TRUE)
  expect_equal(model_comparison$stretch, 2.809, tolerance = 1e-2)
  expect_equal(model_comparison$shift, 0.143, tolerance = 1e-2)
  expect_error(register(gene_data, reference = "Ro19", query = "Col0"))
})
