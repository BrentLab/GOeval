test_that("subset_network runs using total edges", {
  input_file <- file.path("data", "marbach_extract.tsv")
  output_folder <- file.path("data", "marbach_subsets_total")
  name <- "marbach"
  edges <- c(64, 128, 256, 512, 1024, 2048, 4096, 8192)
  subset_network(input_file, output_folder, name, edges)
  expect_equal(length(list.files(file.path("data", "marbach_subsets_total"))), 8)
})

test_that("subset_network runs with num_possible_TFs", {
  input_file <- file.path("data", "marbach_extract.tsv")
  output_folder <- file.path("data", "marbach_subsets_per")
  name <- "marbach"
  edges <- c(16, 32, 64, 128, 256, 512, 1024, 2048)
  num_possible_TFs <- 3
  subset_network(input_file, output_folder, name, edges, num_possible_TFs = num_possible_TFs)
  expect_equal(length(list.files(file.path("data", "marbach_subsets_per"))), 8)
})
