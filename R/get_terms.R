#' Gather ORA results from one folder of ORA summaries
#'
#' @importFrom data.table fread
#'
#' @description
#' Get Over-Representation Analysis (ORA) data for the top n terms from each source node summary file
#'  in a folder output by the `webgestalt_network` function.
#'
#' @param path a folder containing .csv files of the ORA summaries output
#'    by `webgestalt_network`
#' @param n the number of terms from each summary for which to get data
#'
#' @return a data.frame
#'
#' @export
get_terms <- function(path, n) {
  csv_files <- list.files(path = path, pattern = "*.csv", full.names = TRUE)
  if (length(csv_files) > 0) {
    # read the first row of each csv file and combine them into a data frame
    first_rows <- lapply(csv_files, data.table::fread, nrows = n)
    for (i in 1:length(first_rows)) {
      first_rows[[i]]$tfId <- rep(unlist(strsplit(basename(csv_files[[i]]), "_"))[1], each = length(first_rows[[i]]$geneSet))
    }
    combined_df <- do.call(rbind, first_rows)
    # extract TF ID from file name and add to row of corresponding GO term
    # this is now done in `webgestalt_network`
    # combined_df$tfId <- rep(unlist(lapply(csv_files, function(f) {unlist(strsplit(basename(f), '_'))[1]})), each = n)
    return(combined_df)
  } else {
    return(NA)
  }
}
