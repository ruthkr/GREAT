registration_results <- system.file("extdata/brapa_arabidopsis_registration.rds", package = "greatR") |>
  readRDS()

# Summary and visualisation ----
test_that("summarise_registration works", {
  reg_summary <- summarise_registration(registration_results)

  # Expected outputs
  expect_equal(names(reg_summary), c("summary", "registered_genes", "non_registered_genes"))
  expect_equal(length(reg_summary$registered_genes), as.numeric(reg_summary$summary[Result == "Registered genes", Value]))
  expect_equal(length(reg_summary$non_registered_genes), as.numeric(reg_summary$summary[Result == "Non-registered genes", Value]))
})

test_that("plot_registration_results works", {
  gg <- plot_registration_results(registration_results)
  gg_original <- plot_registration_results(registration_results, "original")

  # Expected outputs
  expect_equal(colnames(gg$data), c(colnames(registration_results$data), "gene_facet"))
  expect_equal(nrow(gg$data), nrow(registration_results$data))
  expect_equal(gg$labels$x, "Registered time")
  expect_equal(gg$labels$y, "Scaled expression")
  expect_equal(gg_original$labels$x, "Time point")
  expect_equal(gg_original$labels$y, "Scaled expression")
  expect_no_error(plot_registration_results(registration_results, genes_list = c("BRAA02G018970.3C", "BRAA02G043220.3C")))
  expect_error(plot_registration_results(registration_results, genes_list = 1:2))
})

# Distance and visualisation ----

test_that("calculate_distance works", {
  sample_distance <- calculate_distance(registration_results)

  # Expected outputs
  expect_equal(names(sample_distance), c("registered", "original"))
  expect_equal(unique(sample_distance$registered$timepoint_ref), unique(sample_distance$original$timepoint_ref))
  expect_equal(colnames(sample_distance$registered), colnames(sample_distance$original))
})

test_that("plot_heatmap works", {
  sample_distance <- calculate_distance(registration_results)
  gg <- plot_heatmap(sample_distance, match_timepoints = TRUE)
  gg_original <- plot_heatmap(sample_distance, "original")

  # Expected outputs
  expect_equal(colnames(gg$data), colnames(sample_distance$registered))
  expect_equal(colnames(gg_original$data), colnames(sample_distance$original))
  expect_equal(gg$labels$x, "Col0")
  expect_equal(gg$labels$y, "Ro18")
  expect_equal(gg_original$labels$x, "Col0")
  expect_equal(gg_original$labels$y, "Ro18")
})
