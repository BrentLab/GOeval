#'
#'
#' @importFrom WebGestaltR WebGestaltRBatch
#'
#' @param net
#' @param net_index non-zero integer if net is a permuted version of the real network.
#'  net_index is used in the naming of the output directories.
WG_and_write_results <- function(net, net_index) {
  for (tf in unique_TFs) {
    gene_set = net$V2[net$V1 == tf]
    entrez_gene_set = annotated_ref_genes$entrezgene[annotated_ref_genes$userId %in% gene_set]
    present_terms = annotations_proper_size$geneSet[annotations_proper_size$gene %in% entrez_gene_set]
    if (length(gene_set) > 1 & length(unique(present_terms)) > 1) {
      write.table(gene_set, file=file.path(work_dir,paste0("p",net_index-1),paste0(tf,"_targets.txt")), row.names=FALSE, col.names=FALSE)
    }
  }

  # perform enrichment analysis
  enrich_df <- WebGestaltR::WebGestaltRBatch(
    enrichMethod = METHOD,
    organism = ORGANISM,
    enrichDatabase = DATABASE,
    interestGeneFolder = file.path(work_dir,paste0("p",net_index-1)),
    interestGeneType = GENE_ID,
    referenceGene = annotated_ref_genes$userId,
    referenceGeneType = GENE_ID,
    minNum = 10, # default 10
    maxNum = 500, # default 500
    reportNum = 20, # default 20
    isOutput = GENERATE_REPORT,
    outputDirectory = REPORTS_PATH,
    projectName = paste0(PROJECT_NAME,"_p",net_index-1),
    sigMethod = "top",
    topThr = 16, # top 0.1%
    isParallel = TRUE,
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
        sig_counter = sig_counter + 1
        sig_df['database'] <- rep(DATABASE, nrow(sig_df))
        write.csv(sig_df,file.path(SUMMARIES_PATH,PROJECT_NAME,paste0("p",net_index-1),paste0(tf_method, "_summary.csv")),row.names = FALSE)
      }
    }
  }
  write(paste(sig_counter,"of",length(unique_TFs),"gene sets have at least one significant GO term.\n",length(unique(network$V1)),"gene sets in network."),
        file = file.path(SUMMARIES_PATH,PROJECT_NAME,paste0("p",net_index-1),"README.txt"))
}

## NOTE: running with a previously used PROJECT_NAME will cause overwriting
## NOTE: if using the same PROJECT_NAME for a different network than before,
##         you must first delete the contents of "/webgestalt_work/PROJECT_NAME/"
webgestalt_network <- function(output_directory, network_name, network, reference_set, p) {
  PROJECT_NAME=network_name
  FOLDER_NAME=output_directory
  METHOD="ORA" # ORA | GSEA | NTA
  ORGANISM="scerevisiae" # human: "hsapiens"; yeast: "scerevisiae"
  DATABASE="geneontology_Biological_Process_noRedundant" # see options with listGeneSet()
  GENE_ID="ensembl_gene_id" # see options with listIdType()
  PERMUTATIONS=p # how many randomly permuted networks to create and analyze

  # reports are more in-depth than summaries - advisable to keep reports FALSE if not needed
  REPORTS_PATH="GO_results" # only used if GENERATE_REPORT=TRUE
  SUMMARIES_PATH=file.path("GO_summaries", FOLDER_NAME) # will be created if does not exist
  GENERATE_REPORT=FALSE

  # path must exist even if GENERATE_REPORT=FALSE
  dir.create(REPORTS_PATH, recursive=TRUE)

  for (i in 1:(PERMUTATIONS+1)) {
    dir.create(file.path(SUMMARIES_PATH, PROJECT_NAME, paste0("p",i-1)), recursive=TRUE)
  }

  if (file.exists(NETWORK)) {
    # separate network into files of gene sets
    work_dir = file.path("webgestalt_work", FOLDER_NAME, PROJECT_NAME)
    for (i in 1:(PERMUTATIONS+1)) {
      dir.create(file.path(work_dir, paste0("p",i-1)), recursive=TRUE)
    }

    ref_genes = idMapping(organism = ORGANISM, inputGene = reference_set[1], sourceIdType = GENE_ID, targetIdType = "entrezgene", host = "https://www.webgestalt.org/")
    annotations = loadGeneSet(organism = ORGANISM, enrichDatabase = DATABASE)
    genes_per_term = table(annotations$geneSet$geneSet)
    annotations_proper_size = annotations$geneSet[(genes_per_term[annotations$geneSet$geneSet] >= 10 & genes_per_term[annotations$geneSet$geneSet] <= 500),]
    rm(annotations)
    annotated_ref_genes = ref_genes$mapped[ref_genes$mapped$entrezgene %in% annotations_proper_size$gene,]

    annotated_network = network[(network$V1 %in% annotated_ref_genes$userId & network$V2 %in% annotated_ref_genes$userId),]

    unique_TFs = unique(annotated_network$V1)

    # returns a network with permuted TF and gene labels
    # conserves edge number distributions
    permute <- function(net_mat) {
      permuted_mat <- net_mat
      rownames(permuted_mat) = sample(rownames(permuted_mat))
      colnames(permuted_mat) = sample(colnames(permuted_mat))
      permuted <- data.frame(V1=rep(row.names(permuted_mat), ncol(permuted_mat)), V2=rep(colnames(permuted_mat), each=nrow(permuted_mat)), V3=as.numeric(unlist(permuted_mat)))
      return(permuted[permuted$V3 != 0, 1:2])
    }

    WG_and_write_results(annotated_network, 1)
    ann_network_matrix = xtabs(~ V1 + V2, annotated_network)
    for (i in 2:(PERMUTATIONS+1)) {
      WG_and_write_results(permute(ann_network_matrix), i)
    }

  } else {
    print(paste(NETWORK,"must be an existing file."))
  }
}
