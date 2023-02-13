# brapa_sample_data <- system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "greatR") |>
#   utils::read.csv() |>
#   dplyr::select(
#     gene_id = locus_name,
#     accession,
#     timepoint,
#     replicate = group,
#     expression_value
#   )
#
#
# test_that("get_approximate_stretch works", {
#   approximate_stretch <- get_approximate_stretch(
#     brapa_sample_data,
#     reference = "Ro18",
#     query = "Col0"
#   )
#
#   expect_equal(ceiling(approximate_stretch), 3)
#   expect_equal(floor(approximate_stretch), 2)
#   expect_equal(class(approximate_stretch), "numeric")
# })
#
#
# # test_that("match_names works", {
# #   a <- LETTERS[1:3]
# #   b <- LETTERS[4:5]
# #
# #   results <- match_names(a, lookup = b, error = NULL, name_string = "names", lookup_vec_last = " and ")
# #
# #   expect_error(results, "i.*")
# # })
