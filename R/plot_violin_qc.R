plot_violin_qc <- function(
    seu,
    features = c("percent.mito",
                 "percent.ribo",
                 "percent.globin")
) {
  seu %>%
    Seurat::VlnPlot(pt.size = 0.1, ncol = 3,
                    features = features)
}
