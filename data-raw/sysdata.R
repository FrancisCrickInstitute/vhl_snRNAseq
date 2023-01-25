## code to prepare `sysdata` dataset goes here

# transcript types to get feature percentages for
trancript_types <- readr::read_tsv("data/transcript_types.tsv")

# save to sysdata
usethis::use_data(trancript_types,
                  overwrite = TRUE,
                  internal = TRUE)
