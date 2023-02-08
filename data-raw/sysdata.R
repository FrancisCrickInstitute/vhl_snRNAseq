## code to prepare `sysdata` dataset goes here

# transcript types to get feature percentages for
transcript_types <- readr::read_tsv("data/transcript_types.tsv")

# HIF metagene

# literature markers from PMC6104812
literature_markers <- readr::read_tsv("data/literature_markers.tsv") %>%
  # must fix PVRL4 -> NECTIN
  dplyr::mutate(Gene = dplyr::case_when(Gene == "PVRL4" ~ "NECTIN4",
                                        TRUE ~ Gene)) %>%
  { split(.$Gene, f = as.factor(.$Marker_of)) }

# list of gene modules
gene_modules <- list("tcell" = c("IL7R", "LTB", "TRAC", "CD3D"),
                     "monocyte" = c("CD14", "CST3", "CD68", "CTSS"),
                     "hif" = readr::read_tsv("data/hif_metagene.txt", col_names = F)$X1,
                     pax8 = c("PAX8"))

# save to sysdata
usethis::use_data(transcript_types,
                  literature_markers,
                  gene_modules,
                  overwrite = TRUE,
                  internal = TRUE)
