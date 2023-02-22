registration_results <- system.file("extdata/brapa_arabidopsis_registration.rds", package = "greatR") |>
  readRDS()

# Summary and visualisation ----
test_that("summary_registration works", {
  reg_summary <- summary_registration(registration_results)

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
})
