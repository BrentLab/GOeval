#' Create network subsets with highest scored edges
#'
#' @importFrom dplyr desc
#' @importFrom dplyr arrange
#' @importFrom utils head
#' @importFrom dplyr select
#'
#' @description
#' Given a network with edge scores, `subset_network` sorts the network in descending order of score and
#'  creates subsets using the top given numbers of edges.
#' The size of the subsets can either be determined by the total number of edges (if 'num_possible_TFs' is not specified)
#'  or by the desired average number of target nodes from each possible source node (if 'num_possible_TFs' is set to a number > 0).
#' The latter option multiplies each element of 'edges' by 'num_possible_TFs' to get the number of
#'    edges in each subset.
#'
#' @param input_file a file containing the network to create subsets from. The file should
#'  be tab-separated with three columns: source node, target node, edge score
#' @param output_directory name of the folder in which to place the created network subset files
#' @param name base name for the created network subset files
#' @param edges list of total numbers of edges or average edges per TF to include in each subset
#' @param num_possible_TFs if set to a number > 0, the elements of 'edges' will first be
#'  multiplied by this number to get the number of edges for each subset
#'
#' @details
#' The input file should be tab-separated with three columns: source node (e.g. transcription factor),
#'    target node (e.g. the regulated gene), and edge score.
#' The output files will be tab-separated with two columns: source node, target node.
#'
#' @return NULL
#'
#' @export
subset_network <- function(input_file, output_directory, name, edges, num_possible_TFs = 0) {

  # check existence and format of input_file
  if(!file.exists(input_file)) {
    stop("Network file not found.")
  }

  net_first_line <- tryCatch({
    read.table(input_file, sep = "\t", nrows = 1)
  }, error = function(e) {
    stop("Error reading the network file.")
  })

  net_num_columns <- ncol(net_first_line)
  if (net_num_columns != 2 && net_num_columns != 3) {
    stop("Network file should contain three columns.")
  }

  if(!is.character(name) || grepl(" ", name)) {
    stop("Network name must be a string with no spaces.")
  }

  dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)

  network <- read.table(file = input_file, sep = "\t", header = FALSE)

  if (num_possible_TFs > 0) {
    multiplier <- num_possible_TFs
  } else {
    multiplier <- 1
  }

  # check that there is a third column
  # if not, don't subset
  if (length(colnames(network)) < 3) {
    write.table(network, file.path(output_directory, paste0(name, "_", as.integer(length(network$V1) / multiplier), ".tsv")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    return()
  }

  colnames(network) <- c("Column1", "Column2", "Column3")

  # sort the data frame by the third column in descending order
  sorted_net <- dplyr::arrange(network, dplyr::desc(Column3))

  for (i in edges) {
    if (multiplier * i <= length(sorted_net$Column1)) {
      # get the top rows
      top_rows <- utils::head(sorted_net, multiplier * i)

      # remove the third column
      top_rows <- dplyr::select(top_rows, Column1, Column2)

      # save the top rows as a new .tsv file
      output_file <- paste0(name, "_", i, ".tsv")
      write.table(top_rows, file.path(output_directory, output_file), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
  }
}
