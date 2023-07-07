#' Helper function for `get_metrics`
#'
#' @importFrom WebGestaltR idMapping
#'
#' @description
#' This function calculates the percent of source nodes that are annotated to a GO term
#'  that either regulates or is the same as a GO term for which their target genes are enriched
#'
#' @param terms a data.frame obtained from calling `get_terms` on a folder of the
#'    summaries output by `webgestalt_network`
#' @param organism a string specifying the organism that the data is from, e.g.
#'    "hsapiens" or "scerevisiae"
#' @param gene_id the naming system used for the input genes - see options with WebGestaltR::listIdType()
#'  and see webgestalt.org for examples of each type
#' @param go_ann a data.frame of annotations of genes to GO terms. Obtain with
#'  WebGestaltR::loadGeneSet.
#' @param go_reg a data.frame of the regulatory relationships between GO terms.
#'  Obtain with ontologyIndex::get_ontology.
#'
#' @return a list of source nodes
#'
#' @keywords internal
get_annotation_overlap <- function(terms, organism, gene_id, go_ann, go_reg) {
  tf_ids <- unique(terms$tfId)

  tf_map <- WebGestaltR::idMapping(organism = organism, inputGene = tf_ids, sourceIdType = gene_id, targetIdType = "entrezgene", host = "https://www.webgestalt.org/")

  overlap_list <- list()
  for (tf in tf_ids) {
    # tf_terms = terms[terms$tfId == tf,]
    tf_terms <- terms[(terms$tfId == tf & terms$FDR < 0.05), ]
    if (length(tf_terms$geneSet) > 0) {
      tf_entr <- tf_map$mapped$entrezgene[tf_map$mapped$userId == tf]
      tf_ann <- go_ann[go_ann$gene %in% tf_entr, "geneSet"]
      tf_reg <- go_reg[go_reg$id %in% tf_ann, ]
      reg_list <- c(tf_reg[[2]], tf_reg[[3]], tf_reg[[4]])
      reg_list <- reg_list[reg_list != ""]
      # reg_list = c("")
      if (any(tf_ann %in% tf_terms$geneSet) || any(reg_list %in% tf_terms$geneSet)) {
        overlap_list <- append(overlap_list, tf)
      }
    }
    # sim_df <- data.frame(V1=rep(tf_terms$geneSet, each=length(tf_ann)), V2=rep(tf_ann, length(tf_terms$geneSet)))
    # sim_df <- cbind(sim_df, V3 = apply(sim_df, MARGIN=1, v_goSim))
    # overlap_df <- sim_df[sim_df$V3 > 0.6,]
    # if (nrow(overlap_df) > 0) {
    #  print(overlap_df)
    #  overlap_counter = overlap_counter + 1
    # }
  }

  return(overlap_list)
}

#' Helper function for `get_metrics`
#'
#' @importFrom stats median
#'
#' @description
#' This is a helper function that computes the desired metrics from a data.frame
#'  of the top terms for each source node in a network.
#'
#' @param full_terms a data.frame containing all the summary results for a network.
#'  Use `get_terms` to obtain.
#' @param network_size the number of source nodes in the network
#' @param organism a string specifying the organism that the data is from, e.g.
#'  "hsapiens" or "scerevisiae"
#' @param gene_id the naming system used for the input genes - see options with WebGestaltR::listIdType()
#'  and see webgestalt.org for examples of each type
#' @param get_sum boolean whether to get the 'sum' metric, which is the sum of the negative
#'  log base 10 of the p-value for the top term of each source node minus 'penalty' times the total
#'  number of source nodes.
#' @param get_percent boolean whether to get the 'percent' metric, which is the
#'  percent of source nodes with at least one term with a FDR below the 'fdr_threshold'
#' @param get_mean boolean whether to get the 'mean' metric, which is the mean negative
#'  log base 10 of the p-value for the top term of each source node regardless of significance
#' @param get_median boolean whether to get the 'median' metric, which is the median negative
#'  log base 10 of the p-value for the top term of each source node regardless of significance
#' @param get_annotation_overlap boolean whether to get the 'annotation_overlap' metric,
#'  which is the percent of source nodes that are annotated to at least one of the 16 GO terms for
#'  which their target genes are most enriched
#' @param get_size boolean whether to get the 'size' metric, which is the number of
#'  source nodes in the network subset that have more than one target gene with annotations. This number is
#'  used in the calculation of all other metrics.
#' @param penalty the penalty applied to the 'sum' metric for each source node in the network
#' @param fdr_threshold the FDR threshold for a gene set term to be considered significantly
#'  over-represented for the purposes of calculating the 'percent' metric
#' @param go_ann a data.frame of annotations of genes to GO terms. Obtain with
#'  WebGestaltR::loadGeneSet.
#' @param go_reg a data.frame of the regulatory relationships between GO terms.
#'  Obtain with ontologyIndex::get_ontology.
#'
#' @return a list of metric values
#'
#' @keywords internal
get_network_metrics <- function(full_terms, network_size, organism, gene_id, get_sum, get_percent, get_mean, get_median, get_annotation_overlap, get_size, penalty, fdr_threshold, go_ann = NULL, go_reg = NULL) {
  if (!any(is.na(full_terms))) {
    terms <- full_terms[match(unique(full_terms$tfId), full_terms$tfId), ]
    percent <- 100 * sum(terms$FDR < fdr_threshold) / network_size
    neglogp <- -log10(terms$pValue)
    # cap -logp values due to rounding to 0 for pval < 1E-16
    neglogp <- ifelse(is.finite(neglogp), neglogp, 16)
    #neglogp_zeros <- append(rep(0, network_size - length(terms$geneSet)), neglogp)
    if (get_annotation_overlap) {
      prior_ann <- 100 * length(get_annotation_overlap(full_terms, organism, gene_id, go_ann = go_ann, go_reg = go_reg)) / network_size
    } else {
      prior_ann <- 0
    }
  } else {
    percent <- 0
    neglogp <- rep(0, network_size)
    prior_ann <- 0
  }
  return(c(
    sum(neglogp - penalty), percent, mean(neglogp),
    median(neglogp), prior_ann, network_size
  )[c(get_sum, get_percent, get_mean, get_median, get_annotation_overlap, get_size)])
}

#' Get summary metrics of a network's ORA results
#'
#' @import stringr
#' @importFrom ontologyIndex get_ontology
#' @importFrom WebGestaltR loadGeneSet
#' @importFrom gtools mixedsort
#' @importFrom parallel mclapply
#' @importFrom parallelly availableCores
#'
#' @description
#' `get_metrics` creates a data.frame that contains specified metrics for all network subsets and their permutations
#'  that have a subdirectory in the provided path. It is designed to be run on the output of the `webgestalt_network`
#'  function to prepare summary metrics for plotting with the `plot_metrics` function.
#'  The 'directory' path should contain only directories created by `webgestalt_network`.
#'
#' @param directory a directory containing only the directories of ORA summaries created by
#'  `webgestalt_network` for all networks of interest
#' @param organism a string specifying the organism that the data is from, e.g.
#'  "hsapiens" or "scerevisiae". Only required if get_annotation_overlap = TRUE.
#' @param database the gene set database to search for enrichment - see options with WebGestaltR::listGeneSet().
#'  Must be a Gene Ontology "biological process" database if get_annotation_overlap = TRUE.
#' @param gene_id the naming system used for the input genes - see options with WebGestaltR::listIdType()
#'  and see webgestalt.org for examples of each type. Only required if get_annotation_overlap = TRUE.
#' @param get_sum boolean whether to get the 'sum' metric, which is the sum of the negative
#'  log base 10 of the p-value for the top term of each source node minus 'penalty' times the total
#'  number of source nodes.
#' @param get_percent boolean whether to get the 'percent' metric, which is the
#'  percent of source nodes with at least one term with a FDR below the 'fdr_threshold'
#' @param get_mean boolean whether to get the 'mean' metric, which is the mean negative
#'  log base 10 of the p-value for the top term of each source node regardless of significance
#' @param get_median boolean whether to get the 'median' metric, which is the median negative
#'  log base 10 of the p-value for the top term of each source node regardless of significance
#' @param get_annotation_overlap boolean whether to get the 'annotation_overlap' metric,
#'  which is the percent of source nodes that are annotated to at least one of the 16 GO terms for
#'  which their target genes are most enriched
#' @param get_size boolean whether to get the 'size' metric, which is the number of
#'  source nodes in the network subset that have more than one target gene with annotations. This number is
#'  used in the calculation of all other metrics.
#' @param penalty the penalty applied to the 'sum' metric for each TF in the network
#' @param fdr_threshold the FDR threshold for a gene set term to be considered significantly
#'  over-represented for the purposes of calculating the 'percent' metric
#' @param parallel boolean whether to get the metrics for each network in the directory
#'  in parallel - use with caution, as this has not been adequately tested
#'
#' @return a list of data.frames, each containing the values of one metric.
#'  The columns of a data.frame represent the different subset sizes, and the rows
#'  represent the different network permutations. The first row is from the unpermuted networks.
#'
#' @export
get_metrics <- function(directory, organism = "hsapiens", database = "geneontology_Biological_Process_noRedundant", gene_id = "ensembl_gene_id", get_sum = TRUE, get_percent = FALSE, get_mean = FALSE, get_median = FALSE, get_annotation_overlap = FALSE, get_size = TRUE, penalty = 3, fdr_threshold = 0.05, parallel = FALSE) {
  # GO regulatory relationships only needed if get_annotation_overlap = TRUE
  if (get_annotation_overlap) {
    # Load regulatory relationships between GO terms for the calculation of overlap between TF GO
    # annotations and their targets' enriched GO terms.
    go <- as.data.frame(ontologyIndex::get_ontology(file = "http://purl.obolibrary.org/obo/go/go-basic.obo", extract_tags = "everything"))
    go <- go[go$namespace == 'biological_process',]
    go_reg <- go[c('id', 'regulates', 'negatively_regulates', 'positively_regulates')]
    rm(go)

    # Load gene annotations with WebGestaltR::loadGeneSet.
    suppressWarnings(go_ann <- WebGestaltR::loadGeneSet(organism = organism, enrichDatabase = database, hostName = "https://www.webgestalt.org/")$geneSet)
  } else {
    go_reg <- NULL
    go_ann <- NULL
  }

  networks_list <- list.dirs(directory, full.names = TRUE, recursive = FALSE)
  # filter out the webgestalt_work folder
  networks_list <- networks_list[!grepl("webgestalt_work", networks_list)]

  if (length(networks_list) == 0) {
    stop("Invalid input path: ", directory)
  }

  networks_list <- gtools::mixedsort(networks_list)

  if (parallel) {
    # parallel version
    # works fine when run in debug mode
    # gives "*** recursive gc invocation" when run via the "Test" button
    # but also sometimes it works, so possible race condition although I'm not sure how
    results <- parallel::mclapply(networks_list, function(net) {
      network_data_files <- list.files(net, recursive = TRUE, pattern = "network_data.txt", full.names = TRUE)
      network_sizes <- lapply(network_data_files, function (f) {
        size_line <- grep("valid", readLines(f), value = TRUE)
        return(as.integer(unlist(strsplit(size_line, " "))[1]))
      })
      paths <- list.dirs(net, full.names = TRUE, recursive = FALSE)
      p_terms <- mapply(get_terms, paths, 16, SIMPLIFY = FALSE)
      return(t(mapply(get_network_metrics, p_terms, network_sizes, organism, gene_id, get_sum, get_percent, get_mean, get_median, get_annotation_overlap, get_size, penalty, fdr_threshold, MoreArgs = list(go_ann = go_ann, go_reg = go_reg))))
    }, mc.cores = parallelly::availableCores())
  } else {
    # non-parallel version
    results <- lapply(networks_list, function(net) {
      network_data_files <- list.files(net, recursive = TRUE, pattern = "network_data.txt", full.names = TRUE)
      network_sizes <- lapply(network_data_files, function (f) {
        size_line <- grep("valid", readLines(f), value = TRUE)
        return(as.integer(unlist(strsplit(size_line, " "))[1]))
      })
      paths <- list.dirs(net, full.names = TRUE, recursive = FALSE)
      p_terms <- mapply(get_terms, paths, 16, SIMPLIFY = FALSE)
      return(t(mapply(get_network_metrics, p_terms, network_sizes, organism, gene_id, get_sum, get_percent, get_mean, get_median, get_annotation_overlap, get_size, penalty, fdr_threshold, MoreArgs = list(go_ann = go_ann, go_reg = go_reg))))
    })
  }
  df_list <- list()
  for (i in 1:ncol(results[[1]])) {
    column_list <- lapply(results, "[", , i)
    df_list[[i]] <- data.frame(column_list)
    colnames(df_list[[i]]) <- basename(networks_list)
  }
  return(df_list)
}
