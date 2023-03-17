get_out <- function(out_dir,
                    experiment,
                    genome,
                    sublibrary,
                    parse_analysis_subdir,
                    do_integration,
                    do_timestamp) {

  {
    if (is.null(out_dir))
      # if out_dir not given, use same output structure as in the parse analysis/ directory
      "out/" %>%
      paste(paste(experiment, collapse = "_x_"),
            genome,
            paste(sublibrary, collapse = "_x_"),
            parse_analysis_subdir, sep = "/") %>%
      { if (do_integration) paste0(., "/integrated/") else paste0(., "/unintegrated/") } %>%
      { if (do_timestamp) paste0(., format(Sys.time(), "%Y%m%d_%H%M%S"), "/") else . }
    else
      # if out_dir is given, use out_dir
      out_dir
  } %>% {
    # pre-set file names (TODO: define these once all outputs are finalised)
    purrr::map(list(base = "",
                    cache = "/cache/",
                    infercnv = "/infercnv/",
                    scevan = "/scevan/"),
               function(x)
      paste0(., x, "/")  %>% clean_path())
  }
}
