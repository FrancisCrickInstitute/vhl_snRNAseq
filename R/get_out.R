get_out <- function(out_dir,
                    experiment,
                    genome,
                    sublibrary,
                    parse_analysis_subdir,
                    do_integration,
                    do_timestamp) {
  {
    # if out_dir not given, use same output structure as in the parse analysis/ directory
    if (is.null(out_dir))
      "out/" %>%
      paste(experiment, genome, sublibrary, parse_analysis_subdir, sep = "/") %>%
      { if (do_integration) paste0(., "/integrated/") else paste0(., "/unintegrated/") } %>%
      { if (do_timestamp) paste0(., format(Sys.time(), "%Y%m%d_%H%M%S"), "/") else . }
    # if out_dir is given, use out_dir
    else
      out_dir
  } %>% {
    # pre-set file names (TODO: define these once all outputs are finalised)
    purrr::map(list(base = ""), function(x)
      paste0(., x, "/")  %>% clean_path())
  }
}
