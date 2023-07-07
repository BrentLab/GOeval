
#' Perform the whole GOeval pipeline on one network
#'
#' @description
#' Given one network, `evaluate` will first make subsets with `subset_network`,
#'  then run Over-Representation Analysis (ORA) with `webgestalt_network`,
#'  followed by `get_metrics` and `plot_metrics` to plot selected summary statistics.
#'
#' @param network path to the file containing the network to create subsets from. The file should
#'  be tab-separated with three columns: source node, target node, edge score
#' @param reference_set path to the set of all genes possibly included in the network. Must be a .txt
#'  file containing exactly one column of the genes that could possibly appear in the network.
#' @param output_directory path to the directory in which to store the generated network subsets,
#'  ORA summaries, and plots
#' @param network_name short name for the network - used in file naming so may not contain spaces
#' @param organism a string specifying the organism that the data is from, e.g.
#'    "hsapiens" or "scerevisiae" - see options with WebGestaltR::listOrganism()
#' @param database the gene set database to search for enrichment - see options with WebGestaltR::listGeneSet().
#'  Must be a Gene Ontology "biological process" database if get_annotation_overlap = TRUE.
#' @param gene_id the naming system used for the input genes - see options with WebGestaltR::listIdType()
#'  and see webgestalt.org for examples of each type
#' @param edges list of total numbers of edges or average edges per TF to include in each subset
#' @param num_possible_TFs if set to a number > 0, the elements of 'edges' will first be
#'  multiplied by this number to get the number of edges for each subset
#' @param permutations the number of randomly permuted networks to create and run ORA on
#' @param penalty the penalty applied to the 'sum' metric for each TF in the network
#' @param fdr_threshold the FDR threshold for a gene set term to be considered significantly
#'  over-represented for the purposes of calculating the 'percent' metric
#' @param get_sum boolean whether to get and plot the 'sum' metric, which is the sum of the negative
#'  log base 10 of the p-value for the top term of each source node minus 'penalty' times the total
#'  number of source nodes.
#' @param get_percent boolean whether to get and plot the 'percent' metric, which is the
#'  percent of source nodes with at least one term with a FDR below the 'fdr_threshold'
#' @param get_mean boolean whether to get and plot the 'mean' metric, which is the mean negative
#'  log base 10 of the p-value for the top term of each source node regardless of significance
#' @param get_median boolean whether to get and plot the 'median' metric, which is the median negative
#'  log base 10 of the p-value for the top term of each source node regardless of significance
#' @param get_annotation_overlap boolean whether to get and plot the 'annotation_overlap' metric,
#'  which is the percent of source nodes that are annotated to at least one of the 16 GO terms for
#'  which their target genes are most enriched
#' @param get_size boolean whether to get and plot the 'size' metric, which is the number of
#'  source nodes in the network subset that have more than one target gene with annotations. This number is
#'  used in the calculation of all other metrics.
#' @param plot boolean whether to make plots of the calculated metrics and write them to a pdf
#'
#' @details
#' The input file should be tab-separated with two or three columns: source node (e.g. transcription factor),
#'  target node (e.g. the regulated gene), and, optionally, edge score.
#'
#' @return output of `get_metrics`. Can be used as input to `plot_metrics`.
#'
#' @export
evaluate <- function(network, reference_set, output_directory, network_name, organism = "hsapiens", database = "geneontology_Biological_Process_noRedundant", gene_id = "ensembl_gene_id",
                     edges = c(512, 1024, 2048, 4096, 8192, 16384, 32768, 65536), num_possible_TFs = 0, permutations = 3, penalty = 3, fdr_threshold = 0.05,
                     get_sum = TRUE, get_percent = FALSE, get_mean = FALSE, get_median = FALSE, get_annotation_overlap = FALSE, get_size = TRUE, plot = TRUE) {

  # first package function
  subset_network(network, file.path(output_directory, paste0(network_name, "_subsets")), network_name, edges, num_possible_TFs)

  # for unique names to avoid overwriting
  formatted_time <- format(Sys.time(), "%Y-%m-%d-%H%M")
  summaries_suffix <- paste0("_summaries_", formatted_time)

  # sequential loop over all created subsets
  for (subset in list.files(file.path(output_directory, paste0(network_name, "_subsets")), full.names = TRUE)) {
   # second package function
   webgestalt_network(network_path = subset,
                       reference_set = reference_set,
                       output_directory = file.path(output_directory, paste0(network_name, summaries_suffix)),
                       # this just gets the name of each file minus the extension
                       network_name = strsplit(basename(subset), "[.]")[[1]][1],
                       organism = organism,
                       database = database,
                       gene_id = gene_id,
                       permutations = permutations)
  }

  # third package function
  # mapply returns a list of length 1 containing the output of get_metrics
  metric_dfs_by_net <- mapply(get_metrics, file.path(output_directory, paste0(network_name, summaries_suffix)), MoreArgs=list(organism = organism, database = database, gene_id = gene_id, get_sum = get_sum, get_percent = get_percent, get_mean = get_mean, get_median = get_median, get_annotation_overlap = get_annotation_overlap, get_size = get_size, penalty = penalty, fdr_threshold = fdr_threshold, parallel = FALSE), SIMPLIFY = FALSE)

  if (plot) {
    pdf(file.path(output_directory, paste0(network_name, "_ORA_metrics_plots_", formatted_time, ".pdf")), 7, 5)
    # fourth package function
    plot_data <- plot_metrics(metric_dfs_by_net, network_name, "", perTF = (num_possible_TFs > 0), sum = get_sum, percent = get_percent, mean = get_mean, median = get_median, annotation_overlap = get_annotation_overlap, size = get_size)
    dev.off()
  }

  # return output of get_metrics to easily re-plot with plot_metrics
  return(metric_dfs_by_net)
}
