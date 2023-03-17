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
final_cluster_annotations <- list(
  "221202_A01366_0326_AHHTTWDMXY" = list(
    "1"="tumour",
    "2"="tumour",
    "3"="immune (monocytes)",
    "4"="tumour",
    "5"="endothelial",
    "6"="immune (T cells)"
  )
)

# save to sysdata
usethis::use_data(transcript_types,
                  markers,
                  ditto_colours,
                  umap_void_theme,
                  final_cluster_annotations,
                  overwrite = TRUE,
                  internal = TRUE)
