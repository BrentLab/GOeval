library(ontologyIndex)
suppressMessages(library(WebGestaltR))

test_that("get_metrics and plot_metrics run", {
  #tmp <- system.file("tmp", package = "GOeval")
  #folders <- c(file.path(tmp, "GO_summaries", "np3_npVm_top"), file.path(tmp, "GO_summaries", "marbach_npVm_top"))
  folders <- c(file.path("data", "np3_npVm_top"))

  title_text = "NP3"
  subtitle_text = "subset for overlapping TFs with Marbach et al."

  go <- as.data.frame(get_ontology(file = "http://purl.obolibrary.org/obo/go/go-basic.obo", extract_tags = "everything"))
  go <- go[go$namespace == 'biological_process',]
  go_reg <- go[c('id', 'regulates', 'negatively_regulates', 'positively_regulates')]
  rm(go)

  suppressWarnings(go_ann <- loadGeneSet(organism = "hsapiens", enrichDatabase = "geneontology_Biological_Process_noRedundant", hostName = "https://www.webgestalt.org/")$geneSet)

  #metric_dfs_by_net <- get_metrics(folders, get_percent = TRUE, get_mean = TRUE, get_median = TRUE, get_annotation_overlap = TRUE, go_ann = go_ann, go_reg = go_reg, parallel = FALSE)
  metric_dfs_by_net <- mapply(get_metrics, folders, MoreArgs=list(get_percent = TRUE, get_mean = TRUE, get_median = TRUE, get_annotation_overlap = TRUE, go_ann = go_ann, go_reg = go_reg, parallel = FALSE), SIMPLIFY = FALSE)

  expect_equal(length(metric_dfs_by_net[[1]]), 5)

  #pdf(file.path("data", "test-plot_metrics.pdf"), 7, 5)
  # call plot_metrics
  plot_metrics(metric_dfs_by_net, title_text, subtitle_text)
  #dev.off()

  expect_equal(2 * 2, 4)

  gc()
})
