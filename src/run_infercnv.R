base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
                   "/Volumes/TracerX/",
                   "/camp/project/tracerX/")
setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"))
devtools::load_all()

out_dir <- "out/230210_A01366_0351_AHNHCFDSX5_x_221202_A01366_0326_AHHTTWDMXY/hg38/comb_x_SHE5052A9_S101/all-well/DGE_filtered/unintegrated/"
out <- get_out(out_dir)
dir.create(out$infercnv, showWarnings = F)

# load cds
cds <- readRDS(list.files(out$cache, pattern = "cds\\_singler\\_annotated", full.names = T))

# define lineages
# based on consensus_modules_vs_partitions_heatmap-1.png
final_annotations <- list(
  "1" = "epithelial",
  "2" = "epithelial",
  "3" = "tumour",
  "4" = "endothelial",
  "5" = "tumour",
  "6" = "tumour",
  "7" = "tumour",
  "8" = "myeloid",
  "9" = "tumour",
  "10" = "tumour",
  "11" = "lymphoid",
  "12" = "tumour",
  "13" = "tumour"
)

# get lineages
final_annotations_df <-
  tibble::tibble(partition = names(final_annotations),
                 annotation = unlist(final_annotations)) %>%
  dplyr::mutate(
    lineage = dplyr::case_when(
      annotation %in% c("myeloid", "lymphoid") ~ "immune",
      annotation %in% c("epithelial", "endothelial") ~ "kidney",
      annotation == "tumour" ~ "tumour"
    )
  )
excl_lin <- "immune"
ref_lin <- "kidney"
query_lin <- "tumour"

# encode in the colData
final_annotations_df <-
  dplyr::left_join(
    tibble::as_tibble(SummarizedExperiment::colData(cds)),
    final_annotations_df,
    by = "partition"
  )
SummarizedExperiment::colData(cds) <-
  cbind(
    SummarizedExperiment::colData(cds),
    final_annotations_df
  )

# exclude `excl_lin` cells
cds <- cds[, SummarizedExperiment::colData(cds)$lineage != excl_lin]

# write annotations_file for infercnv
infercnv_annotations_file <- paste0(out$infercnv, "infercnv_annotations.tsv")
SummarizedExperiment::colData(cds)[, "annotation", drop = F] %>%
  write.table(infercnv_annotations_file, col.names = F, quote = F, sep = "\t")

# create infercnv object
infercnv_obj <-
  infercnv::CreateInfercnvObject(
    raw_counts_matrix = cds@assays@data$counts,
    annotations_file = infercnv_annotations_file,
    gene_order_file = "data/gencode/gencode.v43.basic.annotation_clean.bed",
    ref_group_names = c("epithelial", "endothelial")
  )

# perform infercnv operations to reveal cnv signal
infercnv_obj <-
  infercnv::run(
    infercnv_obj,
    cutoff = 1,
    out_dir = out$infercnv,
    cluster_by_groups = T,
    denoise = T,
    HMM = T
  )
