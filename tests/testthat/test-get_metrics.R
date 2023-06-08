test_that("get_metrics runs with annotation overlap", {
  webgestalt_results <- file.path("data", "webgestalt_output", "np3_dtoTFs_summaries")

  metric_dfs <- get_metrics(webgestalt_results, get_annotation_overlap = TRUE, parallel = FALSE)

  expect_equal(length(metric_dfs), 3)

  gc()
})

test_that("get_metrics runs with annotation overlap with only one subset", {
  webgestalt_results <- file.path("data", "webgestalt_output", "EDN_summaries")

  metric_dfs <- get_metrics(webgestalt_results, get_annotation_overlap = TRUE, parallel = FALSE)

  expect_equal(length(metric_dfs), 3)

  gc()
})
