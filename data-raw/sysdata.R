## code to prepare `sysdata` dataset goes here

# transcript types to get feature percentages for
transcript_types <- readr::read_tsv("data/transcript_types.tsv")

# HIF metagene
hif_metagene <- readr::read_tsv("data/hif_metagene.txt", col_names = F)$X1

# literature markers from PMC6104812
literature_markers <- readr::read_tsv("data/literature_markers.tsv") %>%
  # must fix PVRL4 -> NECTIN
  dplyr::mutate(Gene = dplyr::case_when(Gene == "PVRL4" ~ "NECTIN4",
                                        TRUE ~ Gene)) %>%
  { split(.$Gene, f = as.factor(.$Marker_of)) }


# save to sysdata
usethis::use_data(transcript_types,
                  hif_metagene,
                  literature_markers,
                  overwrite = TRUE,
                  internal = TRUE)
