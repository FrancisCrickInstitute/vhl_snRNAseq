# set directory
if (Sys.info()["nodename"]=="Alexs-MacBook-Air-2.local") {
  base_dir <- "/Volumes/TracerX/working/VHL_GERMLINE/tidda/"
} else {
  base_dir <- "/camp/project/tracerX/working/VHL_GERMLINE/tidda/"
}
wkdir <- paste0(base_dir, "/vhl/")
parse_pipeline_dir <- paste0(base_dir, "parse_pipeline")
setwd(wkdir)
library(devtools) ; load_all()

# runs
runs <- paste(
  paste0("( cd ", parse_pipeline_dir, "/analysis/ ;"),
  "find . -mindepth 5 -maxdepth 5 -type d |",
  "grep '.*all-well/.*filtered$' )"
) %>%
  system(intern = T) %>%
  # remove ./ at the beginning of every path
  gsub("^\\.\\/", "", .) %>%
  # convert vectorised bash output to tibble
  { dplyr::tibble(path = .) } %>%
  # convert nested subdirs to variables
  # (each level corresponds to a parameters of the split-pipe run)
  tidyr::separate(
    path,
    into = c("experiment", "genome", "sublibrary", "parse_analysis_subdir"),
    sep = "/",
    remove = F,
    extra = "merge"
  )# %>%
  # add in integration option T/F
  #dplyr::cross_join(dplyr::tibble(do_integration = c(T, F)))
# write to out/
readr::write_tsv(runs, "out/runs.tsv")

base_dir <- get_base_dir()
wkdir <- paste0(base_dir, "/vhl/")
parse_pipeline_dir <- paste0(base_dir, "parse_pipeline")
setwd(wkdir)
runs <- readr::read_tsv("out/runs.tsv") %>%
  dplyr::filter(genome=="hg38", do_integration==F)
for(i in 1:nrow(runs)) {
  run = dplyr::filter(runs, dplyr::row_number() == nrow(runs) + 1 - i)
  cat(run$experiment, run$genome, run$sublibrary, run$parse_analysis_subdir, "\n")
  generate_QC_report(
    experiment = run$experiment,
    parse_pipeline_dir = paste0(base_dir, "/parse_pipeline/"),
    genome = run$genome,
    sublibrary = run$sublibrary,
    parse_analysis_subdir = run$parse_analysis_subdir,
    n_dims = 20,
    out_dir = NULL,
    sample_subset = NULL,
    do_timestamp = F,
    do_integration = run$do_integration)
}

# runs %>% purrr::pwalk(function(...) {
#   run <- tibble::tibble(...)
#   cat(run$experiment, run$genome, run$sublibrary, run$parse_analysis_subdir, "\n")
#   generate_QC_report(
#     experiment = run$experiment,
#     parse_pipeline_dir = "/camp/project/tracerX/working/VHL_GERMLINE/tidda/parse_pipeline/",
#     genome = run$genome,
#     sublibrary = run$sublibrary,
#     parse_analysis_subdir = run$parse_analysis_subdir,
#     n_dims = 20,
#     out_dir = NULL,
#     sample_subset = NULL,
#     do_timestamp = F)
# })
