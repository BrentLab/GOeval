#' plot_metrics
#'
#' Plot the summary metrics from get_metrics across network subset size for one or more networks.
#'
#' @importFrom tidyr gather
#' @importFrom tidyr separate
#' @importFrom dplyr summarise
#' @importFrom dplyr group_by
#' @import ggplot2
#'
#' @param metric_dfs_by_net a list of outputs of get_metrics for each network to be plotted
#' @param title_text text for the title of each plot; generally the names of the networks
#' @param subtitle_text text for the subtitle of each plot
#' @param sum bool whether to plot the 'sum' metric, which is the sum of the negative
#'  log base 10 of the p-value for the top term of each TF minus 3 times the total
#'  number of TFs.
#' @param percent bool whether to plot the 'percent' metric, which is the
#'  percent of TFs with at least one GO term with a FDR < 0.05
#' @param mean bool whether to plot the 'mean' metric, which is the mean negative
#'  log base 10 of the p-value for the top term of each TF
#' @param median bool whether to plot the 'median' metric, which is the median negative
#'  log base 10 of the p-value for the top term of each TF
#' @param annotation_overlap bool whether to plot the 'annotation_overlap' metric,
#'  which is the percent of TFs that are annotated to a GO term for which their
#'  target genes are enriched
#'
#' @export
plot_metrics <- function(metric_dfs_by_net, title_text, subtitle_text, sum = TRUE, percent = FALSE, mean = FALSE, median = FALSE, annotation_overlap = TRUE) {
  to_plot <- c(sum, percent, mean, median, annotation_overlap)
  plot_with_boxes <- FALSE
  for (df_list in metric_dfs_by_net) {
    if (length(df_list[[1]][1,]) == 1) {
      plot_with_boxes <- TRUE
    }
  }

  # merge data frames of like metrics across networks
  metric_dfs <- list()
  for (i in 1:length(metric_dfs_by_net[[1]])) {
    metric_dfs[[i]] <- do.call(cbind, lapply(unname(metric_dfs_by_net), "[[", i))
  }

  label_text <- c("Penalized sum of -logp", "Percent of TFs with significant GO term",
                  "Mean of -logp", "Median of -logp", "Percent of TFs annotated to at least one significant GO term")

  for (m in 1:length(metric_dfs)) {
    gathered = tidyr::gather(metric_dfs[[m]], network, metric)
    gathered <- tidyr::separate(gathered, network, into=c("method", "size"), remove = FALSE)
    gathered$network = factor(gathered$network, levels=unique(gathered$network))

    real = gathered[!duplicated(gathered$network), ]
    real$size <- as.numeric(real$size)
    real$ltype <- "Original network"
    permuted = gathered[duplicated(gathered$network), ]
    permuted$size <- as.numeric(permuted$size)

    plt <- ggplot2::ggplot(real, mapping=ggplot2::aes(x = size, y = metric, color = method, shape = method, linetype = ltype)) + ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::ggtitle(title_text, subtitle_text) +
      ggplot2::guides(linetype = ggplot2::guide_legend(override.aes = list(size = 0.5)),
             color = ggplot2::guide_legend(override.aes = list(size = 0.5))) +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
            panel.background = ggplot2::element_blank(),
            plot.title=ggplot2::element_text(hjust=0.5),
            plot.subtitle=ggplot2::element_text(hjust=0.5)) +
      ggplot2::scale_x_continuous(n.breaks = 10, trans = "log2")

     # Set lower limit of y axis to 0 if graphing a percent or # of TFs
    if (m == 2 || m == 5 || m == 6) {
      plt <- plt + ggplot2::scale_y_continuous(limits = c(0, NA))
    }

    # if any of the networks whose metrics are being plotted has only one subset,
    # plotting with box plots for the permuted networks looks better
    if (plot_with_boxes) {
      plt <- plt + ggplot2::geom_boxplot(permuted, mapping=ggplot2::aes(x=size, y=metric, color = method, group = size))
    } else {
      medians = dplyr::summarise(group_by(permuted, network), median = median(metric))
      medians <- tidyr::separate(medians, network, into=c("method", "size"), remove = FALSE)
      medians$size <- as.numeric(medians$size)
      medians$ltype <- "Randomized network"

      for(net in unique(permuted$method)) {
        plt <- plt + ggplot2::geom_line(medians[medians$method == net,], mapping = ggplot2::aes(x = size, y = median, color = method, linetype = ltype), linewidth = 1.5)
      }
      plt <- plt + ggplot2::labs(x = "Edges per TF", y = label_text[m], color = "Network name", shape = "Network name", linetype = "Data")
    }
    plot(plt)
  }
}
