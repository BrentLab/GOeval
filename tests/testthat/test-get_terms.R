test_that("proper tfId length", {
  tmp <- system.file("tmp", package = "GOeval")
  summaries_path <- file.path(tmp, "GO_summaries", "edn_top", "edn", "p0")
  actual <- get_terms(summaries_path, 1)
  expect_equal(51, length(actual$tfId))
})
