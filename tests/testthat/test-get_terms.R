test_that("proper tfId length", {
  #tmp <- system.file("tmp", package = "GOeval")
  summaries_path <- file.path("data", "np3_npVm_top", "np3_5", "p0")
  actual <- get_terms(summaries_path, 1)
  expect_equal(21, length(actual$tfId))
})
