## code to prepare `sysdata` dataset goes here

# transcript types to get feature percentages for
transcript_types <- readr::read_tsv("data/transcript_types.tsv")

# literature markers
markers <- readRDS("data/markers.rds")

# save to sysdata
usethis::use_data(transcript_types,
                  markers,
                  overwrite = TRUE,
                  internal = TRUE)
