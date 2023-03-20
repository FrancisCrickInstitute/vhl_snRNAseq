library(magrittr)
parse_pipeline_dir <- "../parse_pipeline/"

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
  # (each level corresponds to a parameter of the split-pipe run)
  tidyr::separate(
    path,
    into = c("experiment", "genome", "sublibrary", "parse_analysis_subdir"),
    sep = "/",
    remove = F,
    extra = "merge"
  ) %>%
  # add in integration option T/F
  dplyr::cross_join(dplyr::tibble(do_integration = c(T, F))) %>%
  # if exists, get comb only, else get all SLs
  dplyr::group_by(experiment) %>%
  dplyr::filter(any(sublibrary == "comb") & sublibrary == "comb" |
                  !any(sublibrary == "comb")) %>%
  dplyr::ungroup()

# sample metadata
final_samples <-
  readr::read_tsv(paste0(parse_pipeline_dir,
                         "/expdata/230210_A01366_0351_AHNHCFDSX5/sample_metadata.tsv")) %>%
  # only human samples (removes "^K.*)
  dplyr::filter(grepl("^N.*", sample)) %>%
  dplyr::pull(sample)

# add 8 SLs + 2 SLs runs (all data combined) (vectorised arguments are collapsed with ',')
runs <-
  dplyr::bind_rows(
    runs,
    runs %>%
      dplyr::filter(
        grepl("230210|221202", experiment),
        genome == "hg38") %>%
      dplyr::group_by(genome, do_integration, parse_analysis_subdir) %>%
      dplyr::summarise(dplyr::across(everything(), ~ paste(.x, collapse = ","))) %>%
      dplyr::mutate(sample_subset = paste(final_samples, collapse = ","))
    )

# allocate resources
runs <- runs %>%
  dplyr::mutate(mem = dplyr::case_when(experiment == "221202_A01366_0326_AHHTTWDMXY" ~ 20,
                                       experiment == "230127_A01366_0343_AHGNCVDMXY" ~ 40,
                                       experiment == "230210_A01366_0351_AHNHCFDSX5" ~ 60,
                                       experiment == "230210_A01366_0351_AHNHCFDSX5,221202_A01366_0326_AHHTTWDMXY" ~ 80))

# write to out/
readr::write_tsv(runs, "out/runs.tsv")
