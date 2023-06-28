#' evaluate
#'
#' Perform the whole GOeval pipeline on one network. Given one network, this
#'  function will first make subsets with subset_network, then run ORA with
#'  webgestalt_network, followed by get_metrics and plot_metrics to plot summary
#'  statistics.
#' The input file should be tab-separated with three columns: source node (e.g. transcription factor),
#'  target node (e.g. the regulated gene), and edge score. All genes should be written as Ensembl gene IDs.
#'
#' @param network a file containing the network to create subsets from. The file should
#'  be tab-separated with three columns: source node, target node, edge score
#' @param reference_set path to the set of all genes possibly included in the network. Must be a .txt
#'  file containing exactly one column of the genes that could possibly appear in the network.
#' @param output_directory path to the directory in which to store the generated network subsets,
#'  ORA summaries, and plots
#' @param network_name short name for the network - used in file naming so may not contain spaces
#' @param edges list of total numbers of edges or average edges per TF to include in each subset
#' @param num_possible_TFs if set to a number > 0, the elements of 'edges' will first be
#'  multiplied by this number to get the number of edges for each subset
#' @param permutations the number of randomly permuted networks to create and run ORA
#' @param get_sum bool whether to get the 'sum' metric, which is the sum of the negative
#'  log base 10 of the p-value for the top term of each TF minus 3 times the total
#'  number of TFs.
#' @param get_percent bool whether to get the 'percent' metric, which is the
#'  percent of TFs with at least one GO term with a FDR < 0.05
#' @param get_mean bool whether to get the 'mean' metric, which is the mean negative
#'  log base 10 of the p-value for the top term of each TF
#' @param get_median bool whether to get the 'median' metric, which is the median negative
#'  log base 10 of the p-value for the top term of each TF
#' @param get_annotation_overlap bool whether to get the 'annotation_overlap' metric,
#'  which is the percent of TFs that are annotated to a GO term for which their
#'  target genes are enriched
#' @param get_size bool whether to get the 'size' metric, which is the number of
#'  TFs in the network subset that have more than one target gene with GO annotations
#' @param plot bool whether to make a pdf containing plots of the calculated metrics
#'
#' @export
evaluate <- function(network, reference_set, output_directory, network_name, edges = c(512, 1024, 2048, 4096, 8192, 16384, 32768, 65536), num_possible_TFs = 0,
                     permutations = 3, get_sum = TRUE, get_percent = FALSE, get_mean = FALSE, get_median = FALSE, get_annotation_overlap = FALSE, get_size = TRUE, plot = TRUE) {

  subset_network(network, file.path(output_directory, paste0(network_name, "_subsets")), network_name, edges, num_possible_TFs)

  formatted_time <- format(Sys.time(), "%Y-%m-%d-%H%M")
  summaries_suffix <- paste0("_summaries_", formatted_time)

  for (subset in list.files(file.path(output_directory, paste0(network_name, "_subsets")), full.names = TRUE)) {
   webgestalt_network(network_path = subset,
                       reference_set = reference_set,
                       output_directory = file.path(output_directory, paste0(network_name, summaries_suffix)),
                       # this just gets the name of each file minus the extension
                       network_name = strsplit(basename(subset), "[.]")[[1]][1],
                       permutations = permutations)
  }

  metric_dfs_by_net <- mapply(get_metrics, file.path(output_directory, paste0(network_name, summaries_suffix)), MoreArgs=list(get_sum = get_sum, get_percent = get_percent, get_mean = get_mean, get_median = get_median, get_annotation_overlap = get_annotation_overlap, get_size = get_size, parallel = FALSE), SIMPLIFY = FALSE)

  if (plot) {
    pdf(file.path(output_directory, paste0(network_name, "_ORA_metrics_plots_", formatted_time, ".pdf")), 7, 5)
    plot_data <- plot_metrics(metric_dfs_by_net, network_name, "", perTF = (num_possible_TFs > 0), sum = get_sum, percent = get_percent, mean = get_mean, median = get_median, annotation_overlap = get_annotation_overlap, size = get_size)
    dev.off()
  }

  return(metric_dfs_by_net)
}
