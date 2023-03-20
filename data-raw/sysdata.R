## code to prepare `sysdata` dataset goes here

# transcript types to get feature percentages for
transcript_types <- readr::read_tsv("data/transcript_types.tsv")

# literature markers
markers <- readRDS("data/markers.rds")

# colours
ditto_colours <- list(ggplot2::scale_fill_manual(values = dittoSeq::dittoColors()),
                      ggplot2::scale_colour_manual(values = dittoSeq::dittoColors()))

# umap void theme
umap_void_theme <-
  ggplot2::theme(axis.text = ggplot2::element_blank(),
                 axis.title = ggplot2::element_blank(),
                 axis.ticks = ggplot2::element_blank(),
                 legend.position = "none",
                 plot.margin = ggplot2::unit(c(2,2,2,2), "pt"))

# final cluster annotations
final_annotations_files <-
  list.files("data/final_annotations/", full.names = T)
final_annotations_list <-
  final_annotations_files %>%
  purrr::map(function(fcaf) {
    setNames(readr::read_tsv(fcaf)$annotation,
             readr::read_tsv(fcaf)$group) %>%
      split(readr::read_tsv(fcaf)$lvl)
  }) %>%
  setNames(final_annotations_files %>%
             basename() %>%
             tools::file_path_sans_ext())

# save to sysdata
usethis::use_data(transcript_types,
                  markers,
                  ditto_colours,
                  umap_void_theme,
                  final_annotations_list,
                  overwrite = TRUE,
                  internal = TRUE)
