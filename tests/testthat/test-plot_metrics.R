test_that("get_metrics and plot_metrics run", {
  folders <- c(file.path("data", "np3_dtoTFs_summaries"))

  title_text = "NP3"
  subtitle_text = "subset for overlapping TFs with EDN"

  # call get_metrics
  metric_dfs_by_net <- mapply(get_metrics, folders, MoreArgs=list(get_percent = TRUE, get_mean = TRUE, get_median = TRUE, get_annotation_overlap = TRUE, parallel = FALSE), SIMPLIFY = FALSE)

  expect_equal(length(metric_dfs_by_net[[1]]), 6)

  pdf(file.path("data", "line-test-plot_metrics.pdf"), 7, 5)
  # call plot_metrics
  plot_metrics(metric_dfs_by_net, title_text, subtitle_text, perTF = TRUE, percent = TRUE, mean = TRUE, median = TRUE, annotation_overlap = TRUE)
  dev.off()

  expect_equal(2 * 2, 4)

  gc()
})

test_that("get_metrics and plot_metrics run with single subset network", {
  folders <- c(file.path("data", "np3_dtoTFs_summaries"), file.path("data", "EDN_summaries"))

  title_text = "NP3 vs EDN"
  subtitle_text = "subset for overlapping TFs"

  # call get_metrics
  metric_dfs_by_net <- mapply(get_metrics, folders, MoreArgs=list(get_percent = TRUE, get_mean = TRUE, get_median = TRUE, get_annotation_overlap = TRUE, parallel = FALSE), SIMPLIFY = FALSE)

  expect_equal(length(metric_dfs_by_net[[1]]), 6)

  pdf(file.path("data", "box-test-plot_metrics.pdf"), 7, 5)
  # call plot_metrics
  plot_metrics(metric_dfs_by_net, title_text, subtitle_text, perTF = TRUE, percent = TRUE, mean = TRUE, median = TRUE, annotation_overlap = TRUE)
  dev.off()

  expect_equal(2 * 2, 4)

  gc()
})
