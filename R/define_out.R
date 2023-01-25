define_out <- function(
    experiment,
    genome,
    out_dir,
    do_timestamp
) {
  out <- list(
    base = "",
    args = "arguments_for_analyse_parse.R",
    seu_pre_qc = "seu_pre_qc.rds"
  )
  if (is.null(out_dir)) {
    out <- "out/" %>%
      paste0(experiment, "/", genome, "/") %>%
      { if(do_timestamp) paste0(., format(Sys.time(), "%Y%m%d_%H%M%S"), "/") else . } %>%
      { purrr::map(out, function(x) paste0(., x)) }
  } else {
    out <- out %>% purrr::map(function(x) paste0(out_dir, "/", x))
  }

  dir.create(out$base, showWarnings = F, recursive = T)
  cat("Output will be saved to", out$base, "\n")

  return(out)
}
