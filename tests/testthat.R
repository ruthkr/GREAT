library(testthat)
library(GREAT)

test_check("GREAT", reporter = c("progress", "list", "fail"))
