#' calculate the percent of TFs that are annotated to a GO term for which their
#'  target genes are enriched
#'
#' @importFrom WebGestaltR idMapping
#'
#' @param terms a data.frame obtained from calling get_terms on a folder of the
#'    summaries output by webgestalt_network
#' @param organism a string specifying the organism that the data is from, e.g.
#'    "hsapiens" or "scerevisiae"
#' @param go_ann a data.frame of annotations of genes to GO terms. Obtain with
#'  WebGestaltR::loadGeneSet.
#' @param go_reg a data.frame of the regulatory relationships between GO terms.
#'  Obtain with ontologyIndex::get_ontology.
get_annotation_overlap <- function(terms, organism, go_ann, go_reg) {
  tf_ids <- unique(terms$tfId)

  # assumption of "ensembl_gene_id" is a placeholder
  tf_map <- WebGestaltR::idMapping(organism = organism, inputGene = tf_ids, sourceIdType = "ensembl_gene_id", targetIdType = "entrezgene", host = "https://www.webgestalt.org/")

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

#' helper function to compute the desired metrics from a dataframe of the top
#'  terms for each TF in a single network
#'
#' @importFrom stats median
#'
#' @param full_terms a data.frame containing all the summary results for a network.
#'  Use get_terms to obtain.
#' @param network_size the number of TFs in the network. May be different from
#'  the number of TFs that have summary files.
#' @param organism a string specifying the organism that the data is from, e.g.
#'  "hsapiens" or "scerevisiae". Only specify if get_annotation_overlap = TRUE.
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
#' @param go_ann a data.frame of annotations of genes to GO terms. Obtain with
#'  WebGestaltR::loadGeneSet.
#' @param go_reg a data.frame of the regulatory relationships between GO terms.
#'  Obtain with ontologyIndex::get_ontology.
get_network_metrics <- function(full_terms, network_size, organism, get_sum, get_percent, get_mean, get_median, get_annotation_overlap, go_ann = NULL, go_reg = NULL) {
  # heuristically chosen
  penalty <- 3
  if (!any(is.na(full_terms))) {
    terms <- full_terms[match(unique(full_terms$tfId), full_terms$tfId), ]
    percent <- 100 * sum(terms$FDR < 0.05) / network_size
    neglogp <- -log10(terms$pValue)
    # cap -logp values due to rounding to 0 for pval < 1E-16
    neglogp <- ifelse(is.finite(neglogp), neglogp, 16)
    neglogp_zeros <- append(rep(0, network_size - length(terms$geneSet)), neglogp)
    if (get_annotation_overlap) {
      prior_ann <- 100 * length(get_annotation_overlap(full_terms, organism, go_ann = go_ann, go_reg = go_reg)) / network_size
    } else {
      prior_ann <- 0
    }
  } else {
    percent <- 0
    neglogp_zeros <- rep(0, network_size)
    prior_ann <- 0
  }
  return(c(
    sum(neglogp_zeros - penalty), percent, mean(neglogp_zeros),
    median(neglogp_zeros), prior_ann
  )[c(get_sum, get_percent, get_mean, get_median, get_annotation_overlap)])
}

#' get a data.frame that contains specified metrics for all networks that have a
#'  subdirectory in the provided path including both the real and permuted networks
#'
#' @import stringr
#' @importFrom ontologyIndex get_ontology
#' @importFrom WebGestaltR loadGeneSet
#' @importFrom gtools mixedsort
#' @importFrom parallel mclapply
#' @importFrom reader n.readLines
#' @importFrom parallelly availableCores
#'
#' @param directory a directory containing the webgestalt_network output directory
#'  for each network of interest
#' @param organism a string specifying the organism that the data is from, e.g.
#'    "hsapiens" or "scerevisiae". Only specify if get_annotation_overlap = TRUE.
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
#' @param parallel bool whether to get the metrics for each network in the directory
#'  in parallel or sequentially
#'
#' @export
get_metrics <- function(directory, organism = "hsapiens", get_sum = TRUE, get_percent = FALSE, get_mean = FALSE, get_median = FALSE, get_annotation_overlap = FALSE, parallel = TRUE) {
  # start = Sys.time()

  # GO regulatory relationships only needed if get_annotation_overlap = TRUE
  if (get_annotation_overlap) {
    # Load regulatory relationships between GO terms for the calculation of overlap between TF GO
    # annotations and their targets' enriched GO terms.
    go <- as.data.frame(ontologyIndex::get_ontology(file = "http://purl.obolibrary.org/obo/go/go-basic.obo", extract_tags = "everything"))
    go <- go[go$namespace == 'biological_process',]
    go_reg <- go[c('id', 'regulates', 'negatively_regulates', 'positively_regulates')]
    rm(go)

    # Load gene annotations with WebGestaltR::loadGeneSet.
    suppressWarnings(go_ann <- WebGestaltR::loadGeneSet(organism = "hsapiens", enrichDatabase = "geneontology_Biological_Process_noRedundant", hostName = "https://www.webgestalt.org/")$geneSet)
  } else {
    go_reg <- NULL
    go_ann <- NULL
  }

  networks_list <- list.dirs(directory, full.names = TRUE, recursive = FALSE)
  networks_list <- gtools::mixedsort(networks_list)

  if (parallel) {
    # parallel version
    # works fine when run in debug mode
    # gives "*** recursive gc invocation" when run via the "Test" button
    # but also sometimes it works...
    results <- parallel::mclapply(networks_list, function(net) {
      readme <- reader::n.readLines(file.path(net, "p0", "README.txt"), n = 1)
      network_size <- as.integer(stringr::word(readme, start = 3, end = 3, sep = stringr::fixed(" ")))
      paths <- list.dirs(net, full.names = TRUE, recursive = FALSE)
      p_terms <- mapply(get_terms, paths, 16, SIMPLIFY = FALSE)
      result <- t(mapply(get_network_metrics, p_terms, network_size, organism, get_sum, get_percent, get_mean, get_median, get_annotation_overlap, MoreArgs = list(go_ann = go_ann, go_reg = go_reg)))
      return(result)
    }, mc.cores = parallelly::availableCores())
  } else {
    # non-parallel version
    results <- lapply(networks_list, function(net) {
      readme <- reader::n.readLines(file.path(net, "p0", "README.txt"), n = 1)
      network_size <- as.integer(stringr::word(readme, start = 3, end = 3, sep = stringr::fixed(" ")))
      paths <- list.dirs(net, full.names = TRUE, recursive = FALSE)
      p_terms <- mapply(get_terms, paths, 16, SIMPLIFY = FALSE)
      return(t(mapply(get_network_metrics, p_terms, network_size, organism, get_sum, get_percent, get_mean, get_median, get_annotation_overlap, MoreArgs = list(go_ann = go_ann, go_reg = go_reg))))
    })
  }

  df_list <- list()
  for (i in 1:ncol(results[[1]])) {
    column_list <- lapply(results, "[", , i)
    df_list[[i]] <- data.frame(column_list)
    colnames(df_list[[i]]) <- basename(networks_list)
  }

  # end = Sys.time()
  # print(end - start)
  return(df_list)
}
