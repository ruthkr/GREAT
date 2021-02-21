library(testthat)
library(genereg)

test_check("genereg", reporter = c("progress", "list", "fail"))
