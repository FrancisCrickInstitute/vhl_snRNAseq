plot_scatter_qc <- function(
    seu,
    feature2,
    feature1 = "nCount_RNA"
) {
  seu %>%
    Seurat::FeatureScatter(
      feature1 = feature1,
      feature2 = feature2
    )
}
