library(ontologyIndex)
suppressMessages(library(WebGestaltR))

test_that("get_metrics runs with annotation overlap", {
  #tmp <- system.file("tmp", package = "GOeval")
  webgestalt_results <- file.path("data", "np3_npVm_top")

  # go <- as.data.frame(get_ontology(file = "http://purl.obolibrary.org/obo/go/go-basic.obo", extract_tags = "everything"))
  # go <- go[go$namespace == 'biological_process',]
  # go_reg <- go[c('id', 'regulates', 'negatively_regulates', 'positively_regulates')]
  # rm(go)
  #
  # suppressWarnings(go_ann <- loadGeneSet(organism = "hsapiens", enrichDatabase = "geneontology_Biological_Process_noRedundant", hostName = "https://www.webgestalt.org/")$geneSet)

  metric_dfs <- get_metrics(webgestalt_results, get_annotation_overlap = TRUE, parallel = FALSE)

  expect_equal(length(metric_dfs), 2)

  gc()
})
