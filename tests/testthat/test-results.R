registration_results <- readRDS(system.file("extdata/brapa_arabidopsis_registration.rds", package = "greatR"))

# Summary and visualisation ----
test_that("summary.res_greatR works", {
  reg_summary <- summary(registration_results)

  # Expected outputs
  expect_s3_class(reg_summary, "summary.res_greatR")
  expect_equal(names(reg_summary), c("summary", "registered_genes", "non_registered_genes", "reg_params"))
  expect_equal(length(reg_summary$registered_genes), as.numeric(reg_summary$summary[Result == "Registered genes", Value]))
  expect_equal(length(reg_summary$non_registered_genes), as.numeric(reg_summary$summary[Result == "Non-registered genes", Value]))
})

test_that("plot.summary.res_greatR works", {
  reg_summary <- summary(registration_results)
  gg_all <- plot(reg_summary)
  gg_reg <- plot(reg_summary, type = "registered", type_dist = "density")

  # Expected outputs
  expect_equal(length(gg_all), 4)
  expect_equal(colnames(gg_all[[3]]$data), colnames(reg_summary$reg_params))
  expect_gte(nrow(gg_all[[3]]$data), nrow(gg_reg[[1]]$data))
  expect_equal(gg_all[[3]]$labels$x, "Stretch")
  expect_equal(gg_all[[3]]$labels$y, "Shift")
  expect_equal(gg_all[[1]]$labels$y, "Count")
  expect_equal(gg_all[[4]]$labels$y, "Count")
  expect_equal(gg_reg[[1]]$labels$y, "Density")
  expect_equal(gg_reg[[4]]$labels$y, "Density")
})

test_that("plot.res_greatR works", {
  gg <- plot(registration_results)
  gg_original <- plot(registration_results, "original")

  # Expected outputs
  expect_equal(colnames(gg$data), c(colnames(registration_results$data), "gene_facet"))
  expect_equal(nrow(gg$data), nrow(registration_results$data))
  expect_equal(gg$labels$x, "Registered time")
  expect_equal(gg$labels$y, "Scaled expression")
  expect_equal(gg_original$labels$x, "Time point")
  expect_equal(gg_original$labels$y, "Scaled expression")
  expect_no_error(plot(registration_results, genes_list = c("BRAA02G018970.3C", "BRAA02G043220.3C")))
  expect_error(plot(registration_results, genes_list = 1:2))
})

# Distance and visualisation ----

test_that("calculate_distance works", {
  sample_distance <- calculate_distance(registration_results)

  # Expected outputs
  expect_s3_class(sample_distance, "dist_greatR")
  expect_equal(names(sample_distance), c("result", "original"))
  expect_equal(unique(sample_distance$result$timepoint_ref), unique(sample_distance$original$timepoint_ref))
  expect_equal(colnames(sample_distance$result), colnames(sample_distance$original))
})

test_that("plot.dist_greatR works", {
  sample_distance <- calculate_distance(registration_results)
  gg <- plot(sample_distance, match_timepoints = TRUE)
  gg_original <- plot(sample_distance, "original")

  # Expected outputs
  expect_equal(colnames(gg$data), colnames(sample_distance$result))
  expect_equal(colnames(gg_original$data), colnames(sample_distance$original))
  expect_equal(gg$labels$x, "Col0")
  expect_equal(gg$labels$y, "Ro18")
  expect_equal(gg_original$labels$x, "Col0")
  expect_equal(gg_original$labels$y, "Ro18")
})
