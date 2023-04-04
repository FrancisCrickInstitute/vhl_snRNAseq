# config
base_dir <- ifelse(Sys.info()["nodename"] == "Alexs-MacBook-Air-2.local",
                   "/Volumes/TracerX/",
                   "/camp/project/tracerX/")
setwd(paste0(base_dir,"working/VHL_GERMLINE/tidda/vhl/"))
devtools::load_all()
out_dir <- "out/230210_A01366_0351_AHNHCFDSX5_x_221202_A01366_0326_AHHTTWDMXY/hg38/comb_x_SHE5052A9_S101/all-well/DGE_filtered/unintegrated/"
out <- get_out(out_dir)

# load packages
library(SCPA)
library(msigdbr)
library(magrittr)
library(dplyr)
library(ComplexHeatmap)
library(SummarizedExperiment)
library(circlize)
library(tidyverse)

# functions
fix_pathway_names <- function(pathway_names) {
  pathway_names %>% tolower() %>% gsub("\\_", " ", .) %>% gsub("hallmark ", "", .)
}
plot_my_heatmap <- function(scpa_out_df, col, stat = "qval") {
  scpa_out_df %>%
    select(Pathway, name = matches(col), value = matches(stat)) %>%
    mutate(Pathway = fix_pathway_names(Pathway)) %>%
    pivot_wider() %>%
    column_to_rownames("Pathway") %>%
    as.matrix() %>%
    Heatmap(name = stat,
            border = T,
            show_row_dend = F,
            show_column_dend = T,
            row_names_gp = grid::gpar(fontsize = 8))
}
plot_my_enrichment <- function(p_dat_in, top_n = 5) {
  p_dat <-
    p_dat_in %>%
    mutate(Pathway = Pathway %>% fix_pathway_names(),
           color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                             FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                             FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                             FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))

  ggplot(p_dat, aes(-FC, qval)) +
    geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
    geom_point(cex = 2.6, shape = 21, fill = p_dat$color, stroke = 0.3) +
    ggrepel::geom_text_repel(data = p_dat %>% slice_max(abs(FC), n = top_n), aes(label = Pathway),
                             size = 7, min.segment.length = 0, nudge_x = 10) +
    xlab("Enrichment") +
    ylab("Qval") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          aspect.ratio = 1) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
}

# load cds
cds <-
  readRDS(paste0(out$cache, "cds_annotated.rds"))

# get hallmark pathways
pathways <-
  msigdbr(species = "Homo sapiens", category = "H") %>%
  format_pathways()

# initiate list
scpa_out <- list()

# compare pathways between tumour vs normal
scpa_out[["normal_vs_malignant"]] <-
  compare_sce(cds,
              group1 = "partition_lineage",
              group1_population = c("normal", "malignant"),
              pathways = pathways)

# save
saveRDS(scpa_out[["normal_vs_malignant"]],
        paste0(out$dea, "scpa_normal_vs_malignant.rds"))

# compare pathways between TME-derived normal vs normal-derived normal
# get cell types, collapse lymphoid annots
colData(cds)$cluster_annot <-
  ifelse(grepl("lymphoid", colData(cds)$cluster_annot), "lymphoid",
         colData(cds)$cluster_annot)
cell_types <-
  unique(colData(cds)$cluster_annot)
cell_types <- cell_types[cell_types != "malignant"]

# split tumour and normal
split_cds <-
  list(normal = "normal_renal",
       tumour = c("renal_cyst", "solid", "metastasis")) %>%
  purrr::map(function(lts) {
    cds[, colData(cds)$lesion_type %in% lts]
  })

# run
scpa_out[["normal_tme_vs_normal_normal"]] <- list()
for (ct in cell_types) {
  print(ct)
  normal <- sce_extract(split_cds$normal,
                        meta1 = "cluster_annot",
                        value_meta1 = ct)
  tumour <- sce_extract(split_cds$tumour,
                        meta1 = "cluster_annot",
                        value_meta1 = ct)
  scpa_out[["normal_tme_vs_normal_normal"]][[ct]] <-
    compare_pathways(list(normal, tumour), pathways) %>%
    select(Pathway, qval) %>%
    set_colnames(c("Pathway", paste(ct, "qval", sep = "_")))
}
scpa_out[["normal_tme_vs_normal_normal"]] <-
  scpa_out[["normal_tme_vs_normal_normal"]] %>%
  purrr::reduce(full_join, by = "Pathway") %>%
  set_colnames(gsub(colnames(.), pattern = " ", replacement = "_")) %>%
  select(c("Pathway", grep("_qval", colnames(.)))) %>%
  filter_all(any_vars(. > 2)) %>%
  tibble::column_to_rownames("Pathway")
scpa_out[["normal_tme_vs_normal_normal"]] <-
  scpa_out[["normal_tme_vs_normal_normal"]][rowMeans(scpa_out[["normal_tme_vs_normal_normal"]]) > 0,]
col_hm <- colorRamp2(colors = c("blue", "white", "red"), breaks = c(0, 3, 6))
pdf(paste0(out$dea, "normal_tme_vs_normal_normal_heatmap.pdf"))
Heatmap(scpa_out[["normal_tme_vs_normal_normal"]],
        col = col_hm,
        name = "Qval",
        column_labels = colnames(scpa_out[["normal_tme_vs_normal_normal"]]) %>% gsub("\\_", " ", .) %>% gsub(" qval", "", .),
        row_labels = rownames(scpa_out[["normal_tme_vs_normal_normal"]]) %>% fix_pathway_names(),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        border = T,
        column_km = 3,
        row_km = 2)
dev.off()

list(
  monocle3::plot_cells(
    cds[, colData(cds)$partition_lineage == "normal" &
          colData(cds)$lesion_type != "normal_renal"],
    color_cells_by = "cluster_annot"
  ) + ditto_colours + ggplot2::coord_equal() + ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none"),
  monocle3::plot_cells(
    cds[,colData(cds)$partition_lineage == "normal" &
          colData(cds)$lesion_type == "normal_renal"],
    color_cells_by = "cluster_annot"
  ) + ditto_colours + ggplot2::coord_equal() + ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")
) %>%
  patchwork::wrap_plots(ncol = 2)

# save
saveRDS(scpa_out$normal_tme_vs_normal_normal,
        paste0(out$dea, "normal_tme_vs_normal_normal.rds"))

# normal vs normal at lineage level
# get partition annots
cell_types <-
  unique(colData(cds)$partition_annot)
cell_types <- cell_types[cell_types != "malignant"]
# run
scpa_out[["normal_tme_vs_normal_normal_lineages"]] <- list()
for (ct in cell_types) {
  print(ct)
  normal <- sce_extract(split_cds$normal,
                        meta1 = "partition_annot",
                        value_meta1 = ct)
  tumour <- sce_extract(split_cds$tumour,
                        meta1 = "partition_annot",
                        value_meta1 = ct)
  scpa_out[["normal_tme_vs_normal_normal_lineages"]][[ct]] <-
    compare_pathways(list(normal, tumour), pathways) %>%
    select(Pathway, qval) %>%
    set_colnames(c("Pathway", paste(ct, "qval", sep = "_")))
}
scpa_out[["normal_tme_vs_normal_normal_lineages"]] <-
  scpa_out[["normal_tme_vs_normal_normal_lineages"]] %>%
  purrr::reduce(full_join, by = "Pathway") %>%
  set_colnames(gsub(colnames(.), pattern = " ", replacement = "_")) %>%
  select(c("Pathway", grep("_qval", colnames(.)))) %>%
  filter_all(any_vars(. > 2)) %>%
  tibble::column_to_rownames("Pathway")
# remove pathways that are 0
scpa_out[["normal_tme_vs_normal_normal_lineages"]] <-
  scpa_out[["normal_tme_vs_normal_normal_lineages"]][rowMeans(scpa_out[["normal_tme_vs_normal_normal_lineages"]]) > 0,]
max_q <- ceiling(max(scpa_out[["normal_tme_vs_normal_normal_lineages"]]))
col_hm <- colorRamp2(colors = c("blue", "white", "red"), breaks = c(0, max_q/2, max_q))

pdf(paste0(out$dea, "normal_tme_vs_normal_normal_lineages.pdf"))
Heatmap(as.matrix(scpa_out[["normal_tme_vs_normal_normal_lineages"]]),
        col = col_hm,
        name = "Qval",
        column_labels = colnames(scpa_out[["normal_tme_vs_normal_normal_lineages"]]) %>% gsub("\\_", " ", .) %>% gsub(" qval", "", .),
        row_labels = rownames(scpa_out[["normal_tme_vs_normal_normal_lineages"]]) %>% fix_pathway_names(),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        border = T,
        column_km = 3,
        row_km = 3)
dev.off()

# save rds
saveRDS(scpa_out[["normal_tme_vs_normal_normal_lineages"]],
        paste0(out$dea, "scpa_normal_tme_vs_normal_normal_lineages.rds"))

# comparison of tumour partitions vs normal
malignant_partitions <-
  unique(monocle3::partitions(cds[, SummarizedExperiment::colData(cds)$partition_lineage == "malignant"]))
scpa_out[["normal_vs_malignant_partitions"]] <-
  malignant_partitions %>%
  purrr::map(function(partition) {
    print(partition)
    control <- as.matrix(cds[, SummarizedExperiment::colData(cds)$partition_lineage == "normal"]@assays@data$logcounts)
    test <- as.matrix(cds[, monocle3::partitions(cds) == partition]@assays@data$logcounts)
    compare_pathways(list(control, test), pathways = cancer_pathways)
  }) %>%
  setNames(malignant_partitions) %>%
  dplyr::bind_rows(.id = "partition")

SummarizedExperiment::colData(cds)$normal_or_malig_partition <-
  ifelse(SummarizedExperiment::colData(cds)$partition_lineage == "malignant",
         SummarizedExperiment::colData(cds)$partition,
         ifelse(SummarizedExperiment::colData(cds)$partition_lineage == "normal",
                "normal",
                NA))
monocle3::plot_cells(cds[,!is.na(SummarizedExperiment::colData(cds)$normal_or_malig_partition)], color_cells_by = "normal_or_malig_partition") +
  ditto_colours +
  ggplot2::coord_equal() +
  ggplot2::theme_void()

# save
saveRDS(scpa_out[["normal_vs_malignant_partitions"]],
        paste0(out$dea, "scpa_normal_vs_malignant_partitions.rds"))

# comparison of tumour cells from different lesion types vs normal
malignant_lesion_types <- c("renal_cyst", "metastasis", "solid")
scpa_out[["normal_vs_malignant_lesion_types"]] <-
  malignant_lesion_types %>%
  purrr::map(function(lesion_type) {
    print(lesion_type)
    control <- as.matrix(cds[, SummarizedExperiment::colData(cds)$partition_lineage == "normal"]@assays@data$logcounts)
    test <- as.matrix(cds[, SummarizedExperiment::colData(cds)$lesion_type == lesion_type &
                            SummarizedExperiment::colData(cds)$partition_lineage == "malignant"]@assays@data$logcounts)
    compare_pathways(list(control, test), pathways = pathways)
  }) %>%
  setNames(malignant_lesion_types) %>%
  dplyr::bind_rows(.id = "lesion_type")

# save
saveRDS(scpa_out[["normal_vs_malignant_lesion_types"]],
        paste0(out$dea, "scpa_normal_vs_malignant_lesion_types.rds"))

# get VHL status
cds <-
  readRDS(paste0(out$cache, "cds_infercnv.rds" ))

# comparison normal cells with 3p loss versus without
vhl_statuses <- c("VHL mut", "VHL mut + 3p loss", "3p loss")
scpa_out[["normal_vhl_statuses"]] <-
  vhl_statuses %>%
  purrr::map(function(vhl_stat) {
    print(vhl_stat)
    control <- as.matrix(cds[, SummarizedExperiment::colData(cds)$vhl_status == "WT" &
                               SummarizedExperiment::colData(cds)$partition_lineage == "normal"]@assays@data$logcounts)
    test <- as.matrix(cds[, SummarizedExperiment::colData(cds)$vhl_status == vhl_stat &
                            SummarizedExperiment::colData(cds)$partition_lineage == "normal"]@assays@data$logcounts)
    compare_pathways(list(control, test), pathways = pathways)
  }) %>%
  setNames(vhl_statuses) %>%
  dplyr::bind_rows(.id = "vhl_status")

# save
saveRDS(scpa_out[["normal_vhl_statuses"]],
        paste0(out$dea, "scpa_normal_vhl_statuses.rds"))

# vhl loss in epithelial cells
control <- as.matrix(cds[, SummarizedExperiment::colData(cds)$vhl_status == "WT" &
                           SummarizedExperiment::colData(cds)$partition_annot == "epithelial"]@assays@data$logcounts)
chr3_loss_subset <-
  SummarizedExperiment::colData(cds)$has_loss_chr3 == 1 &
  SummarizedExperiment::colData(cds)$partition_annot == "epithelial"
chr3_loss_subset[is.na(chr3_loss_subset)] <- FALSE
test <- as.matrix(cds[, chr3_loss_subset]@assays@data$logcounts)
scpa_out[["normal_chr3_loss"]] <-
    compare_pathways(list(control, test), pathways = pathways)

# find top markers
top_chr3_loss_markers <-
  monocle3::top_markers(cds[c(rownames(control), rownames(test)),],
                        group_cells_by = "has_loss_chr3") %>%
  dplyr::arrange(cell_group, desc(pseudo_R2)) %>%
  dplyr::group_by(cell_group) %>%
  dplyr::top_n(n=5, pseudo_R2)

# comparison of met and cyst vs solid
malignant_nonsolid_lesion_types <- c("renal_cyst", "metastasis")
scpa_out[["malignant_solid_vs_malignant_nonsolid"]] <-
  malignant_nonsolid_lesion_types %>%
  purrr::map(function(lesion_type) {
    print(lesion_type)
    control <- as.matrix(cds[, SummarizedExperiment::colData(cds)$lesion_type == "solid" &
                               SummarizedExperiment::colData(cds)$partition_lineage == "malignant"]@assays@data$logcounts)
    test <- as.matrix(cds[, SummarizedExperiment::colData(cds)$lesion_type == lesion_type &
                            SummarizedExperiment::colData(cds)$partition_lineage == "malignant"]@assays@data$logcounts)
    compare_pathways(list(control, test), pathways = pathways)
  }) %>%
  setNames(malignant_nonsolid_lesion_types) %>%
  dplyr::bind_rows(.id = "lesion_type")

# save
saveRDS(scpa_out[["malignant_solid_vs_malignant_nonsolid"]],
        paste0(out$dea, "scpa_malignant_solid_vs_malignant_nonsolid.rds"))

# read
scpa_files <- list.files(out$dea, pattern = "scpa", full.names = T)
scpa_out <-
  scpa_files %>%
  lapply(readRDS) %>%
  setNames(scpa_files %>% basename %>% tools::file_path_sans_ext() %>% gsub("scpa\\_", "", .))

# get clinical and genotype annotations
clin_data <-
  readr::read_tsv(paste0(out$base, "clinical_and_genotype_annotations.tsv"))
growth_rates <-
  readr::read_tsv(paste0(out$base, "growth_rates.tsv"))

plot_my_heatmap(scpa_out$normal_vhl_statuses, "vhl_status")
plot_my_heatmap(scpa_out$normal_vs_malignant_partitions, "partition")
plot_my_heatmap(scpa_out$normal_vs_malignant_lesion_types, "lesion_type")
plot_my_heatmap(scpa_out$malignant_solid_vs_malignant_nonsolid, "lesion_type")
plot_my_enrichment(scpa_out$malignant_solid_vs_malignant_nonsolid)
plot_my_enrichment(scpa_out$normal_chr3_loss,
                   top_n = 6)
plot_my_enrichment(scpa_out$malignant_solid_vs_malignant_nonsolid %>%
                     filter(lesion_type == "metastasis"),
                   top_n = 5)
plot_my_enrichment(scpa_out$malignant_solid_vs_malignant_nonsolid %>%
                     filter(lesion_type == "renal_cyst"),
                   top_n = 5)
plot_my_enrichment(scpa_out$normal_vs_malignant_lesion_types %>%
                     filter(lesion_type == "solid"),
                   top_n = 5)
plot_my_enrichment(scpa_out$normal_vs_malignant_lesion_types %>%
                     filter(lesion_type == "renal_cyst"),
                   top_n = 5)
plot_my_enrichment(scpa_out$normal_vs_malignant,
                   top_n = 5)
dev.off()

