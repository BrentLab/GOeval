#' plot the summary metrics from get_metrics across network subset size for one or more networks
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
  #print("got here 0")
  to_plot <- c(sum, percent, mean, median, annotation_overlap)
  plot_with_boxes <- FALSE
  for (df_list in metric_dfs_by_net) {
    if (length(df_list[[1]][1,]) == 1) {
      plot_with_boxes <- TRUE
    }
  }
  #print("got here 1")
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
    permuted = gathered[duplicated(gathered$network), ]
    permuted$size <- as.numeric(permuted$size)
    #print("got here 2")
    # if any of the networks whose metrics are being plotted has only one subset,
    # plotting with box plots for the permuted networks looks better
    if (plot_with_boxes) {
      box_plot <- ggplot2::ggplot(real, mapping=ggplot2::aes(x = size, y = metric, color = method, group = method)) + ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point(real, mapping=ggplot2::aes(x = size, y = metric, color = method, shape = method), size = 2.5) +
        ggplot2::ggtitle(title_text, subtitle_text) +
        ggplot2::labs(x = "Edges per TF", y = label_text[m]) +
        ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 0.5)),
               color = ggplot2::guide_legend(override.aes = list(size = 0.5))) +
        ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
              panel.background = ggplot2::element_blank(),
              plot.title=ggplot2::element_text(hjust=0.5),
              plot.subtitle=ggplot2::element_text(hjust=0.5)) +
        ggplot2::scale_x_continuous(n.breaks = 10, trans = "log2")
      #print("got here 3")
      box_plot <- box_plot + ggplot2::geom_boxplot(permuted, mapping=ggplot2::aes(x=size, y=metric, color = method, group = size))
      plot(box_plot)
      #print("got here 4")
    } else {
      medians = dplyr::summarise(group_by(permuted, network), median = median(metric))
      medians <- tidyr::separate(medians, network, into=c("method", "size"), remove = FALSE)
      medians$size <- as.numeric(medians$size)

      line_plot <- ggplot2::ggplot(real, mapping=ggplot2::aes(x = size, y = metric, color = method, group = method)) + ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point(real, mapping=ggplot2::aes(x = size, y = metric, color = method, shape = method), size = 2.5) +
        ggplot2::ggtitle(title_text, subtitle_text) +
        ggplot2::labs(x = "Edges", y = label_text[m]) +
        ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 0.5)),
               color = ggplot2::guide_legend(override.aes = list(size = 0.5))) +
        ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
              panel.background = ggplot2::element_blank(),
              plot.title=ggplot2::element_text(hjust=0.5),
              plot.subtitle=ggplot2::element_text(hjust=0.5)) +
        ggplot2::scale_x_continuous(n.breaks = 10, trans = "log2")

      for(net in unique(permuted$method)) {
        line_plot <- line_plot + ggplot2::geom_line(medians[medians$method == net,], mapping = ggplot2::aes(x = size, y = median, color = method, group = method), linewidth = 1.5, linetype="dotted")
      }
      plot(line_plot)
    }
  }


  # if (plot_with_boxes) {
  #   for (m in 1:length(metric_dfs)) {
  #     if(to_plot[[m]]) {
  #       gathered = tidyr::gather(metric_dfs[[m]], network, metric)
  #       gathered <- tidyr::separate(gathered, network, into=c("method", "size"), remove = FALSE)
  #       gathered$network = factor(gathered$network, levels=unique(gathered$network))
  #
  #       box_plot <- ggplot2::ggplot(gathered, ggplot2::aes(x = network, y = metric, color = method)) +
  #         ggplot2::labs(x = "Edges per TF", y = label_text[m]) +
  #         ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 0.5)),
  #                         color = ggplot2::guide_legend(override.aes = list(size = 0.5))) +
  #         ggplot2::theme(panel.background = ggplot2::element_blank()) +
  #         ggplot2::scale_x_discrete(labels = gathered$size[c(TRUE, diff(as.numeric(gathered$size)) != 0)])
  #
  #       plot(box_plot + ggplot2::geom_boxplot() +
  #              ggplot2::stat_summary(ggplot2::aes(x = network, y = metric), fun=first, geom = "point", color = "red"))
  #     }
  #   }
  # } else {
  #   for (m in 1:length(metric_dfs)) {
  #     if(to_plot[[m]]) {
  #       gathered = tidyr::gather(metric_dfs[[m]], network, metric)
  #       gathered <- tidyr::separate(gathered, network, into=c("method", "size"), remove = FALSE)
  #       gathered$network = factor(gathered$network, levels=unique(gathered$network))
  #
  #       gathered$size = factor(gathered$size, levels=unique(gathered$size))
  #       real = gathered[!duplicated(gathered$network), ]
  #       permuted = gathered[duplicated(gathered$network), ]
  #
  #       line_plot <- ggplot2::ggplot() + ggplot2::geom_line(real, mapping=ggplot2::aes(x = size, y = metric, color = method, group = method), linewidth = 1) +
  #         ggplot2::labs(x = "Edges per TF", y = label_text[m]) +
  #         ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 0.5)),
  #                         color = ggplot2::guide_legend(override.aes = list(size = 0.5))) +
  #         ggplot2::theme(panel.background = ggplot2::element_blank())
  #
  #       medians = dplyr::summarise(dplyr::group_by(permuted, network), median = median(metric))
  #       medians <- tidyr::separate(medians, network, into=c("method", "size"), remove = FALSE)
  #
  #       plot(line_plot + ggplot2::geom_boxplot(permuted, mapping=ggplot2::aes(x=size, y=metric, color = method)) +
  #              ggplot2::geom_point(real, mapping=ggplot2::aes(x = size, y = metric, color = method, shape = method), size = 2.5))
  #     }
  #   }
  # }
}
