library(ontologyIndex)
suppressMessages(library(WebGestaltR))

test_that("get_metrics runs with annotation overlap", {
  webgestalt_results <- file.path("data", "np3_npVm_top")

  metric_dfs <- get_metrics(webgestalt_results, get_annotation_overlap = TRUE, parallel = FALSE)

  expect_equal(length(metric_dfs), 2)

  gc()
})
