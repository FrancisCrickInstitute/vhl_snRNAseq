base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
                   "/Volumes/TracerX/",
                   "/camp/project/tracerX/")
setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"))
devtools::load_all()

out_dir <- "out/230210_A01366_0351_AHNHCFDSX5_x_221202_A01366_0326_AHHTTWDMXY/hg38/comb_x_SHE5052A9_S101/all-well/DGE_filtered/unintegrated/"

# final markers table
final_markers <-
  readr::read_tsv("../renal_sc/data/processed/markers/markers.tsv")

final_markers %>%
  tidyr::separate_longer_delim(study, ", ") %>%
  dplyr::arrange(population) %>%
  dplyr::group_by(population) %>%
  dplyr::summarise(
    n = dplyr::n(),
    gene = paste(unique(gene), collapse = ", "),
    study = paste(unique(study), collapse = ", ")
  ) %>%
  dplyr::arrange(desc(n)) %>%
  readr::write_tsv(paste0(out_dir, "final_markers.tsv"))
