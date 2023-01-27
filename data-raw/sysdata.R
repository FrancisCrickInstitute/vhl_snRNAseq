## code to prepare `sysdata` dataset goes here

# transcript types to get feature percentages for
transcript_types <- readr::read_tsv("data/transcript_types.tsv")

# HIF metagene
hif_metagene <- readr::read_tsv("data/hif_metagene.txt", col_names = F)$X1

# save to sysdata
usethis::use_data(transcript_types,
                  hif_metagene,
                  overwrite = TRUE,
                  internal = TRUE)
