## code to prepare `sysdata` dataset goes here

# transcript types to get feature percentages for
transcript_types <- readr::read_tsv("data/transcript_types.tsv")

# literature markers from PMC6104812
kidney_lineage_markers <- readr::read_tsv("data/kidney_lineage_markers.tsv") %>%
  # must fix PVRL4 -> NECTIN
  dplyr::mutate(Gene = dplyr::case_when(Gene == "PVRL4" ~ "NECTIN4",
                                        TRUE ~ Gene)) %>%
  { split(.$Gene, f = as.factor(.$Marker_of)) }

# list of gene modules
# markers from Dropbox/Parse_pilot/initial_QC_and_analysis.Rmd (curated by Annika Fendler)
celltype_markers <- list(hif = readr::read_tsv("data/hif_metagene.txt", col_names = F)$X1,
                     monocyte = c("CD14", "CST3", "CD68", "CTSS"),
                     endothelial = c("PECAM1", "CLDN5", "ERG", "CDH5", "CD34"),
                     macrophage = c("AIF1", "CD68", "LST1", "IFITM2", "PADI2", "ITGAM"),
                     dendritic = c("THBD", "CLEC9A", "CLEC4C"),
                     b = c("CD79A", "MS4A1", "LINC00926"),
                     collagen_biosynthesis = c("COL1A2", "PDGFA", "PDGFRA", "SMA"),
                     ccrcc = c("BAP1", "PBRM1", "SETD2", "VHL"),
                     fibroblast = c("COL1A1", "COL1A2", "COL5A1", "LOXL1", "LUM", "FBLN1", "FBLN2"),
                     t = c("CD8A", "CD3D", "CD4", "CD3E", "NCAM1", "PTPRC"),
                     nk = c("NCAM1", "NCR1"),
                     pax8 = c("PAX8"))

# save to sysdata
usethis::use_data(transcript_types,
                  kidney_lineage_markers,
                  celltype_markers,
                  overwrite = TRUE,
                  internal = TRUE)
