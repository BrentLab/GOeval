#' webgestalt_network
#'
#' Given a network or network subset, webgestalt_network will generate ORA results
#' for this network and a given number of permutations.
#' The input network file should be a tab-separated .tsv where the first column is the
#' source nodes (TFs) and the second column is the target nodes (regulated genes).
#' All genes must be written as their Ensembl gene ID.
#' NOTE: if using the same output_directory and network_name for a different network
#'  than before, you must first delete the contents of "/webgestalt_work/network_name/"
#'
#' @importFrom WebGestaltR idMapping
#' @importFrom WebGestaltR loadGeneSet
#' @importFrom WebGestaltR WebGestaltRBatch
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom utils write.csv
#' @importFrom stats xtabs
#' @importFrom parallelly availableCores
#'
#' @param network_path path to the network or network subset to evaluate. Must be a tab-separated
#'  file where the first column is the source nodes (TFs) and the second column is the target nodes (regulated genes).
#' @param reference_set path to the set of all genes possibly included in the network. Must be a .txt
#'  file containing exactly one column of the genes that could possibly appear in the network.
#' @param output_directory path to the folder in which output from all networks subset from the same original network should be stored
#' @param network_name the name of the folder to store the results within the output_directory.
#'  It should be in the same format as the files output by subset_network: "\{name\}_\{# of edges\}" i.e. "example_8".
#' @param organism human: "hsapiens"; yeast: "scerevisiae"
#' @param database the gene set database to search for enrichment - see options with WebGestaltR::listGeneSet()
#' @param permutations the number of randomly permuted networks to create and run ORA
#'
#' @export
webgestalt_network <- function(network_path, reference_set, output_directory, network_name, organism = "hsapiens", database = "geneontology_Biological_Process_noRedundant", permutations = 10) {
  METHOD="ORA" # ORA | GSEA | NTA
  GENE_ID="ensembl_gene_id" # see options with listIdType() - set to ensembl_gene_id for now
  # reports are more in-depth than summaries - advisable to keep reports FALSE if not needed
  REPORTS_PATH=output_directory#"GO_results" # only used if GENERATE_REPORT=TRUE
  GENERATE_REPORT=FALSE

  # path must exist even if GENERATE_REPORT=FALSE
  dir.create(REPORTS_PATH, showWarnings = FALSE, recursive=TRUE)

  for (i in 1:(permutations+1)) {
    dir.create(file.path(output_directory, network_name, paste0("p",i-1)), showWarnings = FALSE, recursive=TRUE)
  }

  if (file.exists(network_path)) {
    # create paths for the intermediary gene set files
    work_dir = file.path(output_directory, "webgestalt_work", network_name)
    for (i in 1:(permutations+1)) {
      dir.create(file.path(work_dir, paste0("p",i-1)), showWarnings = FALSE, recursive=TRUE)
    }

    # finds the subset of genes in the reference set that are annotated to valid GO terms
    ref_genes = WebGestaltR::idMapping(organism = organism, inputGene = read.table(reference_set)$V1, sourceIdType = GENE_ID, targetIdType = "entrezgene", host = "https://www.webgestalt.org/")
    annotations = suppressWarnings(WebGestaltR::loadGeneSet(organism = organism, enrichDatabase = database))
    genes_per_term = table(annotations$geneSet$geneSet)
    annotations_proper_size = annotations$geneSet[(genes_per_term[annotations$geneSet$geneSet] >= 10 & genes_per_term[annotations$geneSet$geneSet] <= 500),]
    rm(annotations)
    annotated_ref_genes = ref_genes$mapped[ref_genes$mapped$entrezgene %in% annotations_proper_size$gene,]

    network = read.table(file=network_path, sep='\t', header=FALSE)
    annotated_network = network[(network$V1 %in% annotated_ref_genes$userId & network$V2 %in% annotated_ref_genes$userId),]

    # returns a network with permuted TF and gene labels
    # conserves edge number distributions
    permute <- function(net_mat) {
      permuted_mat <- net_mat
      rownames(permuted_mat) = sample(rownames(permuted_mat))
      colnames(permuted_mat) = sample(colnames(permuted_mat))
      permuted <- data.frame(V1=rep(row.names(permuted_mat), ncol(permuted_mat)), V2=rep(colnames(permuted_mat), each=nrow(permuted_mat)), V3=as.numeric(unlist(permuted_mat)))
      return(permuted[permuted$V3 != 0, 1:2])
    }

    unique_TFs <- unique(annotated_network$V1)
    ann_network_matrix = xtabs(~ V1 + V2, annotated_network)

    # separates each network subset into gene set files, calls the WebGestaltRBatch
    # function on that folder, and saves the results
    for (i in 1:(permutations+1)) {
      if (i == 1) {
        net <- annotated_network
      } else {
        net <- permute(ann_network_matrix)
      }
      valid_sets_counter = 0
      for (tf in unique_TFs) {
        gene_set = net$V2[net$V1 == tf]
        entrez_gene_set = annotated_ref_genes$entrezgene[annotated_ref_genes$userId %in% gene_set]
        present_terms = annotations_proper_size$geneSet[annotations_proper_size$gene %in% entrez_gene_set]
        if (length(gene_set) > 1 & length(unique(present_terms)) > 1) {
          write.table(gene_set, file=file.path(work_dir,paste0("p",i-1),paste0(tf,"_targets.txt")), row.names=FALSE, col.names=FALSE)
          valid_sets_counter = valid_sets_counter + 1
        }
      }

      # perform ORA
      enrich_df <- WebGestaltR::WebGestaltRBatch(
        enrichMethod = METHOD,
        organism = organism,
        enrichDatabase = database,
        interestGeneFolder = file.path(work_dir,paste0("p",i-1)),
        interestGeneType = GENE_ID,
        referenceGene = annotated_ref_genes$userId,
        referenceGeneType = GENE_ID,
        minNum = 10, # default 10
        maxNum = 500, # default 500
        reportNum = 20, # default 20
        isOutput = GENERATE_REPORT,
        outputDirectory = REPORTS_PATH,
        projectName = paste0(network_name,"_p",i-1),
        sigMethod = "top",
        topThr = 16, # top 0.1%
        isParallel = TRUE,
        nThreads = ifelse(nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_", "")) && Sys.getenv("_R_CHECK_LIMIT_CORES_", "") == "TRUE", 2, parallelly::availableCores()),
        hostName = "https://www.webgestalt.org/"
      )

      sig_counter = 0
      for (set in enrich_df) {
        # get name from input
        file_name = unlist(strsplit(set[[1]], '/'))[length(unlist(strsplit(set[[1]], '/')))]
        no_ext = unlist(strsplit(file_name, '[.]'))[1]
        tf_method = paste0(unlist(strsplit(no_ext, '_'))[1], '_', METHOD)

        # save summary as .csv files
        if (!is.null(set[[2]])) {
          df <- data.frame(set[[2]])
          sig_df <- subset(df, select = -c(link))
          if (nrow(sig_df) > 0) {
            if (sig_df$FDR[1] < 0.05) {
              sig_counter = sig_counter + 1
            }
            sig_df['database'] <- rep(database, nrow(sig_df))
            write.csv(sig_df,file.path(output_directory,network_name,paste0("p",i-1),paste0(tf_method, "_summary.csv")),row.names = FALSE)
          }
        }
      }
      write(paste0(sig_counter," significant\n",valid_sets_counter," valid\n",length(unique(network$V1))," total\n",length(network$V2)," edges"),
            file = file.path(output_directory,network_name,paste0("p",i-1),"network_data.txt"))
    }
  } else {
    print(paste(network_path,"must be an existing file."))
  }
}
