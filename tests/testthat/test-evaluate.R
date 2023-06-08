test_that("evaluate runs on 10,000 edge network", {
  evaluate(file.path("data", "marbach_extract.tsv"), file.path("data", "marbach_NP_subset_universe.txt"), "data",
           "mar", permutations = 1, get_annotation_overlap = TRUE)
  expect_equal(2 * 2, 4)
})
