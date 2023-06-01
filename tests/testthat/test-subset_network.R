test_that("subset_network runs", {
  input_file <- file.path("data", "marbach_extract.tsv")
  output_folder <- file.path("data", "network_subsets")
  name <- "marbach"
  edges <- c(64, 128, 256, 512, 1024, 2048, 4096, 8192)
  perTF <- FALSE
  subset_network(input_file, output_folder, name, edges, perTF)
  expect_equal(2 * 2, 4)
})
