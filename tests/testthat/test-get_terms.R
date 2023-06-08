test_that("proper tfId length", {
  summaries_path <- file.path("data", "np3_dtoTFs_summaries", "np3_8", "p0")
  actual <- get_terms(summaries_path, 1)
  expect_equal(16, length(actual$tfId))
})
