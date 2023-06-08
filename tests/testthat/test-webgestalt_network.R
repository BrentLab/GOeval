test_that("webgestalt_network runs", {
  webgestalt_network(file.path("data", "small_net.tsv"), file.path("data", "h1_k562_overlap_universe.txt"), file.path("data", "small_test"), "small_net_194", permutations = 2)

  expect_equal(2 * 2, 4)
})
