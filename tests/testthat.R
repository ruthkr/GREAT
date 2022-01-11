library(testthat)
library(greatR)

test_check("greatR", reporter = c("progress", "list", "fail"))
