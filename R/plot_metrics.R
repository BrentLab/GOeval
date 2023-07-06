#' Plot the summary metrics from `get_metrics`
#'
#' @importFrom tidyr gather
#' @importFrom tidyr separate_wider_regex
#' @importFrom dplyr summarise
#' @importFrom dplyr group_by
#' @import ggplot2
#'
#' @description
#' `plot_metrics` creates plots of all the summary metrics calculated by `get_metrics`.
#' By making the 'metric_dfs_by_net' argument a list of multiple outputs from `get_metrics`,
#'  you can plot the metrics from multiple networks on the same graphs.
#' The x-axis is based on the network subset sizes.
#'
#' @param metric_dfs_by_net a list of outputs from `get_metrics`
#' @param title_text text for the title of every plot; generally the names of the networks
#' @param subtitle_text text for the subtitle of every plot
#' @param perTF boolean whether the network subset sizes were specified as average target
#'  genes per TF; changes the x-axis label
#' @param sum boolean whether to plot the 'sum' metric, which is the sum of the negative
#'  log base 10 of the p-value for the top term of each source node minus the penalty times the total
#'  number of source nodes.
#' @param percent boolean whether to plot the 'percent' metric, which is the
#'  percent of source nodes with at least one term with a FDR below a threshold
#' @param mean boolean whether to plot the 'mean' metric, which is the mean negative
#'  log base 10 of the p-value for the top term of each source node regardless of significance
#' @param median boolean whether to plot the 'median' metric, which is the median negative
#'  log base 10 of the p-value for the top term of each source node regardless of significance
#' @param annotation_overlap boolean whether to plot the 'annotation_overlap' metric,
#'  which is the percent of source nodes that are annotated to at least one of the 16 GO terms for
#'  which their target genes are most enriched
#' @param size boolean whether to plot the 'size' metric, which is the number of
#'  source nodes in the network subset that have more than one target gene with annotations
#'
#' @return a list of data.frames each containing values of one metric.
#'  The column names denote the network and subset. The first row of each data.frame
#'  is from the real networks; the rest are from permuted networks.
#'
#' @export
plot_metrics <- function(metric_dfs_by_net, title_text, subtitle_text = "", perTF, sum = TRUE, percent = FALSE, mean = FALSE, median = FALSE, annotation_overlap = FALSE, size = TRUE) {
  # determine if any of the input networks have only one size
  # if yes, metrics for permuted networks will be shown with box plots instead of line plots
  plot_with_boxes <- FALSE
  if (inherits(metric_dfs_by_net[[1]], "list")) {
    for (df_list in metric_dfs_by_net) {
      if (length(df_list[[1]][1,]) == 1) {
        plot_with_boxes <- TRUE
      }
    }
  } else if (length(metric_dfs_by_net[[1]][1,]) == 1) {
    plot_with_boxes <- TRUE
    metric_dfs_by_net <- list(metric_dfs_by_net)
  } else {
    metric_dfs_by_net <- list(metric_dfs_by_net)
  }

  # merge data frames of like metrics across networks
  metric_dfs <- list()
  for (i in 1:length(metric_dfs_by_net[[1]])) {
    metric_dfs[[i]] <- do.call(cbind, lapply(unname(metric_dfs_by_net), "[[", i))
  }

  # assemble y-axis labels
  to_plot <- c(sum, percent, mean, median, annotation_overlap, size)
  label_text <- c("Sum of (top -log(p-value) - 3) across TFs", "Percent of TFs with a significant GO term (FDR < 0.05)",
                  "Mean of top -log(p-value)", "Median of top -log(p-value)", "Percent of TFs annotated to a significant GO term (FDR < 0.05)",
                  "Number of TFs with > 1 annotated target gene")[to_plot]

  # this won't catch when correct amount but wrong selection of metrics
  if (length(label_text) != length(metric_dfs)) {
    stop("The selected metrics must be the same as those calculated by get_metrics")
  }

  for (m in 1:length(metric_dfs)) {
    # gather dfs into format usable by ggplot2
    gathered = tidyr::gather(metric_dfs[[m]], network, metric)
    gathered <- tidyr::separate_wider_regex(gathered, network, c(method = ".*", "_", size = ".*"), cols_remove = FALSE)
    gathered$network = factor(gathered$network, levels=unique(gathered$network))

    # split metrics based on whether they came from real or permuted data
    real = gathered[!duplicated(gathered$network), ]
    real$size <- as.numeric(real$size)
    real$ltype <- "Original network"
    permuted = gathered[duplicated(gathered$network), ]
    permuted$size <- as.numeric(permuted$size)

    if (plot_with_boxes) {
      plt <- ggplot2::ggplot(real, mapping=ggplot2::aes(x = size, y = metric, color = method, shape = method))
    } else {
      plt <- ggplot2::ggplot(real, mapping=ggplot2::aes(x = size, y = metric, color = method, shape = method, linetype = ltype))
    }

    plt <- plt + ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::ggtitle(title_text, subtitle_text) +
      ggplot2::guides(linetype = ggplot2::guide_legend(override.aes = list(size = 0.5)),
             color = ggplot2::guide_legend(override.aes = list(size = 0.5))) +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
            panel.background = ggplot2::element_blank(),
            plot.title=ggplot2::element_text(hjust=0.5),
            plot.subtitle=ggplot2::element_text(hjust=0.5)) +
      ggplot2::scale_x_continuous(n.breaks = 10, trans = "log2")

     # set lower limit of y axis to 0 if graphing a percent
    if (grepl("Percent", label_text[m])) {
      plt <- plt + ggplot2::scale_y_continuous(limits = c(0, NA))
    }

    # set lower limit of y axis to 0 and ensure integer values if graphing number of source nodes
    # breaks function from https://stackoverflow.com/questions/15622001/how-to-display-only-integer-values-on-an-axis-using-ggplot2
    if (grepl("Number", label_text[m])) {
      plt <- plt + ggplot2::scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))), limits = c(0, NA))
    }

    # if any of the networks whose metrics are being plotted has only one size,
    # plotting with box plots for the permuted networks looks better
    if (plot_with_boxes) {
      plt <- plt + ggplot2::geom_boxplot(permuted, mapping=ggplot2::aes(x=size, y=metric, color = method, group = size)) +
        ggplot2::labs(x = ifelse(perTF, "Edges per TF", "Edges"), y = label_text[m], color = "Network name", shape = "Network name")
    } else {
      medians = dplyr::summarise(group_by(permuted, network), median = median(metric))
      medians <- tidyr::separate_wider_regex(medians, network, c(method = ".*", "_", size = ".*"), cols_remove = FALSE)
      medians$size <- as.numeric(medians$size)
      medians$ltype <- "Randomized network"

      for(net in unique(permuted$method)) {
        plt <- plt + ggplot2::geom_line(medians[medians$method == net,], mapping = ggplot2::aes(x = size, y = median, color = method, linetype = ltype), linewidth = 1.5)
      }
      plt <- plt + ggplot2::labs(x = ifelse(perTF, "Edges per TF", "Edges"), y = label_text[m], color = "Network name", shape = "Network name", linetype = "Data")
    }
    plot(plt)
  }
  # each element of metric_dfs is a DataFrame containing values of one metric
  # the column names denote the network and subset
  # the first row is from the real data; the rest are from the permuted data
  return(metric_dfs)
}
