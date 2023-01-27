boxplot_top_genes <- function(seu, n_genes = 20) {
  # matrix of raw counts
  cts <- seu %>%
    Seurat::GetAssayData(assay = "RNA", slot = "counts") %>%
    as.matrix()

  # get percentage/cell
  cts <- t(cts) / colSums(cts) * 100
  medians <- apply(cts, 2, median)

  # get top n genes
  most_expressed <- order(medians, decreasing = T)[n_genes:1]
  most_exp_matrix <- as.matrix((cts[, most_expressed]))

  # prepare for plotting
  most_exp_df <- stack(as.data.frame(most_exp_matrix))
  colnames(most_exp_df) <- c("perc_total", "gene")

  # boxplot with ggplot2
  boxplot <-
    ggplot2::ggplot(most_exp_df, ggplot2::aes(x = gene, y = perc_total)) +
    ggplot2::geom_boxplot() +
    ggplot2::coord_flip()

  return(boxplot)
}
