get_out <- function(out_dir = NULL,
                    out_subdir = NULL,
                    experiment = NULL,
                    genome = NULL,
                    sublibrary = NULL,
                    parse_analysis_subdir = NULL,
                    do_integration = F,
                    do_timestamp = F) {

  if (is.null(out_dir)) {
    collapsed_names <-
      tibble::tibble(experiment = experiment, sublibrary = sublibrary) %>%
      dplyr::arrange(dplyr::desc(experiment)) %>%
      dplyr::summarise(dplyr::across(everything(), ~ paste(.x, collapse = "_x_")))
  }

  {
    if (is.null(out_dir))
      # if out_dir not given, use same output structure as in the parse analysis/ directory
      # standardise combined experiments - sort and combine experiment and sublibrary by descending exp order
      "out/" %>%
      paste(collapsed_names$experiment,
            genome,
            collapsed_names$sublibrary,
            parse_analysis_subdir,
            sep = "/") %>%
      { if (do_integration) paste0(., "/integrated/") else paste0(., "/unintegrated/") }
    else
      # if out_dir is given, use out_dir
      out_dir
  } %>% {
    if (!is.null(out_subdir))
      paste0(., "/", out_subdir, "/")
    else
      .
  } %>% {
    # timestamp
    if (do_timestamp)
      paste0(., format(Sys.time(), "%Y%m%d_%H%M%S"), "/")
    else
      .
  } %>% {
    # pre-set file names (TODO: define these once all outputs are finalised)
    purrr::map(list(base = "",
                    cache = "/cache/",
                    infercnv = "/infercnv/",
                    scevan = "/scevan/",
                    dea = "/dea/"),
               function(x)
      paste0(., x, "/")  %>% clean_path())
  }

}
