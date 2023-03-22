#metric_dfs_by_net <- mapply(get_metric_dfs, FOLDERS, MoreArgs=list(go_ann, go_reg), SIMPLIFY = FALSE)

plot_metrics <- function(metric_dfs_by_net, figure_name) {
  metric_dfs <- list()
  for (i in 1:length(metric_dfs_by_net[[1]])) {
    metric_dfs[[i]] <- do.call(cbind, lapply(unname(metric_dfs_by_net), "[[", i))
  }

  label_text <- c("Penalized sum of -logp", "Percent of TFs with significant GO term",
                  "Mean of -logp", "Median of -logp", "Percent of TFs annotated to at least one significant GO term")

  pdf(file.path("figures", figure_name), 7, 5)
  for (m in 1:length(metric_dfs)) {
    gathered = gather(metric_dfs[[m]], network, metric)
    gathered <- separate(gathered, network, into=c("method", "size"), remove = FALSE)
    gathered$network = factor(gathered$network, levels=unique(gathered$network))

    box_plot <- ggplot(gathered, aes(x = network, y = metric, color = method)) +
      labs(x = "Edges per TF", y = label_text[m]) +
      guides(shape = guide_legend(override.aes = list(size = 0.5)),
             color = guide_legend(override.aes = list(size = 0.5))) +
      theme(panel.background = element_blank()) +
      scale_x_discrete(labels = gathered$size[c(TRUE, diff(as.numeric(gathered$size)) != 0)])

    plot(box_plot + geom_boxplot() +
           stat_summary(aes(x = network, y = metric), fun=first, geom = "point", color = "red"))

    gathered$size = factor(gathered$size, levels=unique(gathered$size))
    real = gathered[!duplicated(gathered$network), ]
    permuted = gathered[duplicated(gathered$network), ]

    line_plot <- ggplot() + geom_line(real, mapping=aes(x = size, y = metric, color = method, group = method), linewidth = 1) +
      labs(x = "Edges per TF", y = label_text[m]) +
      guides(shape = guide_legend(override.aes = list(size = 0.5)),
             color = guide_legend(override.aes = list(size = 0.5))) +
      theme(panel.background = element_blank())

    medians = summarise(group_by(permuted, network), median = median(metric))
    medians <- separate(medians, network, into=c("method", "size"), remove = FALSE)

    plot(line_plot + geom_boxplot(permuted, mapping=aes(x=size, y=metric, color = method)) +
           geom_point(real, mapping=aes(x = size, y = metric, color = method, shape = method), size = 2.5))
  }
  dev.off()
}
