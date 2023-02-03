define_out <- function(
    experiment,
    genome,
    sublibrary,
    parse_analysis_subdir,
    out_dir,
    do_timestamp,
    do_integration
) {
  out <- list(
    base = "",
    args = "arguments_for_analyse_parse.R",
    seu_pre_qc = "seu_pre_qc.rds",
    seu_post_filtering = "seu_post_filtering.rds",
    seu_transformed = "seu_tranformed.rds",
    seu_split_transformed = "seu_split_transformed.rds",
    seu_umap = "seu_umap.rds",
    seu_celldex_annot = "seu_celldex_annot.rds"
  )
  if (is.null(out_dir)) {
    out <- "out/" %>%
      paste(experiment, genome, sublibrary, parse_analysis_subdir, sep = "/") %>%
      { if(do_integration) paste0(., "/integrated/") else . } %>%
      { if(do_timestamp) paste0(., format(Sys.time(), "%Y%m%d_%H%M%S"), "/") else . } %>%
      { purrr::map(out, function(x) paste0(., x)) }
  } else {
    out <- out %>% purrr::map(function(x) paste0(out_dir, "/", x))
  }

  dir.create(out$base, showWarnings = F, recursive = T)
  cat("Output will be saved to", out$base, "\n")

  return(out)
}
