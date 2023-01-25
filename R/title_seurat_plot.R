title_seurat_plot <- function(
    plotlist,
    title
) {
  # adds title to seurat-generated plots
  # must input plotlist (return.plotlist = T or combine = F)
  p <- cowplot::plot_grid(plotlist = plotlist)
  title <- cowplot::ggdraw() + cowplot::draw_label(title, fontface = 'bold')
  out <- cowplot::plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1))
  return(out)
}
