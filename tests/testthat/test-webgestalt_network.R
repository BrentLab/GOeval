test_that("webgestalt_network runs", {
  tmp <- system.file("tmp", package = "GOeval")
  webgestalt_network(file.path(tmp, "small_net.tsv"), file.path(tmp, "h1_k562_overlap_universe.txt"), file.path(tmp, "small_test"), "small_net_3", permutations = 2)

  expect_equal(2 * 2, 4)
})
