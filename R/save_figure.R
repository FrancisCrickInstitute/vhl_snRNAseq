save_figure <- function(plots, name, out_dir, width, height, res, type = "pdf"){
  if (type == "png") {
    png(paste0(out_dir, name, ".", type),
        width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(out_dir, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}
