## code to prepare `sysdata` dataset goes here

# transcript types to get feature percentages for
transcript_types <- readr::read_tsv("data/transcript_types.tsv")

# literature markers
markers <- readRDS("data/markers.rds")

# hif metagene
hif_metagene <- list(hif = readr::read_tsv("data/hif_metagene.txt", col_names = F)$X1)

# save to sysdata
usethis::use_data(transcript_types,
                  markers,
                  hif_metagene,
                  overwrite = TRUE,
                  internal = TRUE)
