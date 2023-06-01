#' subset_network
#'
#' From a network with edge scores, sorts the network in descending order of score, and
#'  creates subsets using the top given numbers of edges.
#'  The input file should be tab-separated with three columns: source node (e.g. transcription factor),
#'    target node (e.g. the regulated gene), and edge score.
#'  The output files will be tab-separated with two columns: source node, target node.
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr n_distinct
#' @importFrom dplyr select
#'
#' @param input_file a file containing the network to create subsets from. The file should
#'  be tab-separated with three columns: source node, target node, edge score
#' @param output_directory name of the folder in which to place the created network subset files
#' @param name base name for the created network subset files
#' @param edges list of numbers of edges or average edges per TF to include in each subset
#' @param perTF boolean TRUE if the elements of "edges" represent average edges per TF
#'  in the subsets and FALSE if they represent the total number of edges
#'
#' @export
subset_network <- function(input_file, output_directory, name, edges, perTF) {
  network = read.table(file=input_file, sep='\t', header=FALSE)

  colnames(network) <- c("Column1", "Column2", "Column3")

  # Sort the data frame by the third column in descending order
  sorted_net <- dplyr::arrange(network, desc(Column3))

  for (i in edges) {
    # Get the top rows for each file
    top_rows <- head(sorted_net, ifelse(perTF, dplyr::n_distinct(sorted_net$Column1)*i, i))

    # Remove the third column
    top_rows <- dplyr::select(top_rows, Column1, Column2)

    # Save the top rows as a new .tsv file
    output_file <- paste0(name, "_", i, ".tsv")
    write.table(top_rows, file.path(output_directory, output_file), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
}
