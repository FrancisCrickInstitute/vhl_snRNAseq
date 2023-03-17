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
final_cluster_annotations_files <-
  list.files("data/", pattern = "final_cluster_annotations", full.names = T)
final_cluster_annotations_list <-
  final_cluster_annotations_files %>%
  purrr::map(function(fcaf) {
    setNames(readr::read_tsv(fcaf, col_names = F)$X2,
             readr::read_tsv(fcaf, col_names = F)$X1)
  }) %>%
  setNames(final_cluster_annotations_files %>%
             basename() %>%
             tools::file_path_sans_ext() %>%
             gsub("final_cluster_annotations_", "", .))

# gene order file for inferCNV
system(
  paste(
    "mkdir data/gencode ; cd data/gencode",
    "wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz",
    "gunzip gencode.v43.basic.annotation.gtf.gz",
    "cat gencode.v43.basic.annotation.gtf | awk 'OFS=\"\t\" {if ($3==\"gene\") {print $1,$4-1,$5,$14}}' | tr -d '\";' > gencode.v43.basic.annotation.bed",
    sep = ";"
  )
)
readr::read_tsv("data/gencode/gencode.v43.basic.annotation.bed",
                col_names = c("chr", "start", "end", "gene"),
                show_col_types = F) %>%
  # remove weird chromosomes
  dplyr::filter(!greplany(c("fix", "alt", "random", "chrUn", "chrM"), chr)) %>%
  # collapse to full range of each gene
  dplyr::group_by(gene) %>%
  dplyr::summarize(chr = unique(chr),
                 start = min(start),
                 end = max(end)) %>%
  # remove duplicate genes
  dplyr::filter(dplyr::n() == 1) %>%
  readr::write_tsv("data/gencode/gencode.v43.basic.annotation_clean.bed", col_names = F)

# save to sysdata
usethis::use_data(transcript_types,
                  markers,
                  ditto_colours,
                  umap_void_theme,
                  final_cluster_annotations_list,
                  overwrite = TRUE,
                  internal = TRUE)
