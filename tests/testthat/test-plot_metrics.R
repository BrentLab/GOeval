library(ontologyIndex)
suppressMessages(library(WebGestaltR))

test_that("get_metrics and plot_metrics run", {
  folders <- c(file.path("data", "np3_npVm_top")) # file.path("..", "..", "tmp", "GO_summaries", "marbach_npVm_top"),

  title_text = "NP3"
  subtitle_text = "subset for overlapping TFs with Marbach et al."

  # call get_metrics
  metric_dfs_by_net <- mapply(get_metrics, folders, MoreArgs=list(get_percent = TRUE, get_mean = TRUE, get_median = TRUE, get_annotation_overlap = TRUE, parallel = FALSE), SIMPLIFY = FALSE)

  expect_equal(length(metric_dfs_by_net[[1]]), 5)

  pdf(file.path("data", "test-plot_metrics.pdf"), 7, 5)
  # call plot_metrics
  plot_metrics(metric_dfs_by_net, title_text, subtitle_text)
  dev.off()

  expect_equal(2 * 2, 4)

  gc()
})
